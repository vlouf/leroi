#include <bom/io/configuration.h>
#include <bom/io/cf.h>
#include <bom/io/nc.h>
#include <bom/io/odim.h>
#include <bom/radar/beam_propagation.h>
#include <bom/array2.h>
#include <bom/ellipsoid.h>
#include <bom/grid_coordinates.h>
#include <bom/grid_transform.h>
#include <bom/map_projection.h>
#include <bom/trace.h>
#include <getopt.h>

#include <algorithm>
#include <filesystem>
#include <fstream>

using namespace bom;

constexpr auto example_config =
R"(# example layered-flow config

# domain projection
proj4 "+proj=aea +lat_1=-32.2 +lat_2=-35.2 +lon_0=151.209 +lat_0=-33.7008 +a=6378137 +b=6356752.31414 +units=m"

# grid size
size "401 401"

# top left coordinates
left_top "-100250 100250"

# grid resolution
cell_delta "500 -500"

# horizontal grid units
units m

# altitude of lowest layer (m)
altitude_base 100.0

# altitude step between layers (m)
altitude_step 200.0

# number of layers
layer_count 40

# radar moment to generate CAPPIs from
moment DBZH

# whether to output the cappis as well as flow fields
output_cappis true

# whether to output the flow magnitude and angle fields
output_polar true

# maximum distance from CAPPI altitude to use reflectivities
max_alt_dist 20000

# exponent for inverse distance weighting when interpolating between vertical levels (2 is a good default)
idw_pwr 2.0
)";

// Argument parser
constexpr auto try_again = "try --help for usage instructions\n";
constexpr auto usage_string =
R"(LeROI

usage:
  leroi [options] vol.h5 out.nc

available options:
  -h, --help
      Show this message and exit
  -g, --generate
      Output a sample configuration file and exit
)";
constexpr auto short_options = "hg";
constexpr struct option long_options[] =
{
    {"help", no_argument, 0, 'h'},
    {"generate", no_argument, 0, 'g'},
    {0, 0, 0}
};


#if 1
constexpr float nodata = std::numeric_limits<float>::quiet_NaN();
constexpr float undetect = -32.0f;
#else
constexpr float nodata = -32.0f;
constexpr float undetect = -32.0f;
#endif


struct bin_info
{
  float slant_range;
  float ground_range;
  float altitude;
};


struct sweep
{
  radar::beam_propagation beam;
  array1<bin_info>        bins; // @ bin centers
  array1<angle>           rays; // @ ray centers
  array2f                 data;
};


struct volume
{
  latlonalt     location;
  vector<sweep> sweeps;
};


auto find_ground_range_bin(array1<bin_info> const& bins, float target) -> size_t
{
    auto ipos = std::upper_bound(bins.begin(), bins.end(), target, [](auto& l, auto& r) { return l < r.ground_range; });
    if (ipos == bins.end())
        return bins.size();
    if (ipos == bins.begin())
        return 0;
    if (std::fabs(ipos->ground_range - target) < std::fabs((ipos - 1)->ground_range - target))
        return ipos - bins.begin();
    return ipos - bins.begin() - 1;
}


auto find_ray(array1<angle> const& rays, angle target) -> size_t{
    // HACK - not handling 0/360 well
    auto ipos = std::upper_bound(rays.begin(), rays.end(), target);
    if (ipos == rays.end())
        return rays.size() - 1;
    if (ipos == rays.begin())
        return 0;
    if (abs(*ipos - target) < abs(*(ipos - 1) - target))
        return ipos - rays.begin();
    return ipos - rays.begin() - 1;
}


auto generate_cappi(
    volume const& vol,
    array2<latlon> const& latlons,
    float max_alt_diff,
    float idw_pwr,
    float altitude
) -> array2f{
    auto cappi = array2f{latlons.extents()};

    for (size_t y = 0; y < cappi.extents().y; ++y){
        for (size_t x = 0; x < cappi.extents().x; ++x){
            auto ll = latlons[y][x];
            auto br = wgs84.latlon_to_bearing_range(vol.location, ll);
            br.first = br.first.normalize();

            int lwr_scan = -1, upr_scan = -1;
            int lwr_bin, upr_bin;
            float lwr_val, upr_val;
            float lwr_dist, upr_dist;

            for (size_t iscan = 0; iscan < vol.sweeps.size(); ++iscan){
                auto& scan = vol.sweeps[iscan];

                auto ibin = find_ground_range_bin(scan.bins, br.second);
                if (ibin >= scan.bins.size())
                    continue;

                auto alt_dist = scan.bins[ibin].altitude - altitude;
                if (std::fabs(alt_dist) > max_alt_diff)
                    continue;

                auto iray = find_ray(scan.rays, br.first);
                auto val = scan.data[iray][ibin];
                if (std::isnan(val))
                    continue;

                if (alt_dist <= 0){
                    if (lwr_scan == -1 || scan.bins[ibin].altitude > vol.sweeps[lwr_scan].bins[lwr_bin].altitude){
                        lwr_scan = iscan;
                        lwr_bin = ibin;
                        lwr_val = val;
                        lwr_dist = -alt_dist;
                    }
                }
                else{
                    if (upr_scan == -1 || scan.bins[ibin].altitude < vol.sweeps[upr_scan].bins[upr_bin].altitude){
                        upr_scan = iscan;
                        upr_bin = ibin;
                        upr_val = val;
                        upr_dist = alt_dist;
                    }
                }
            }

            if (lwr_scan != -1 && upr_scan != -1){
                if (lwr_dist > 0.0f){
                    auto idw_lwr = 1.0 / std::pow(lwr_dist, idw_pwr);
                    auto idw_upr = 1.0 / std::pow(upr_dist, idw_pwr);
                    auto norm = idw_lwr + idw_upr;
                    cappi[y][x] = lwr_val * (idw_lwr / norm) + upr_val * (idw_upr / norm);
                }
                else
                    cappi[y][x] = lwr_val;
            }
            else if (lwr_scan != -1)
                cappi[y][x] = lwr_val;
            else if (upr_scan != -1)
                cappi[y][x] = upr_val;
            else
                cappi[y][x] = nodata;
        }
    }
    return cappi;
}


auto init_altitudes(io::configuration const& config) -> array1f{
    auto alts = array1f{config["layer_count"]};
    auto base = float(config["altitude_base"]);
    auto step = float(config["altitude_step"]);
    for (auto i = 0; i < alts.size(); ++i)
        alts[i] = base + step * i;

    return alts;
}


auto read_volume(std::filesystem::path const& path, string moment) -> volume
{
    auto vol_odim = io::odim::polar_volume{path, io_mode::read_only};
    auto vol = volume{};

    vol.location.lat = vol_odim.latitude() * 1_deg;
    vol.location.lon = vol_odim.longitude() * 1_deg;
    vol.location.alt = vol_odim.height();

    for (size_t iscan = 0; iscan < vol_odim.scan_count(); ++iscan){
        auto scan_odim = vol_odim.scan_open(iscan);
        auto scan = sweep{};

        scan.beam = radar::beam_propagation{vol.location.alt, scan_odim.elevation_angle() * 1_deg};

        scan.bins.resize(scan_odim.bin_count());
        auto range_scale = scan_odim.range_scale();
        auto range_start = scan_odim.range_start() * 1000 + range_scale * 0.5;
        for (size_t i = 0; i < scan.bins.size(); ++i){
            scan.bins[i].slant_range = range_start + i * range_scale;
            std::tie(scan.bins[i].ground_range, scan.bins[i].altitude) = scan.beam.ground_range_altitude(scan.bins[i].slant_range);
        }

        scan.rays.resize(scan_odim.ray_count());
        auto ray_scale = 360_deg / scan.rays.size();
        auto ray_start = scan_odim.ray_start() * 1_deg + ray_scale * 0.5;
        for (size_t i = 0; i < scan.rays.size(); ++i)
            scan.rays[i] = ray_start + i * ray_scale;

        for (size_t idata = 0; idata < scan_odim.data_count(); ++idata){
            auto data_odim = scan_odim.data_open(idata);
            if (data_odim.quantity() != moment)
                continue;

            scan.data.resize(vec2z{scan_odim.bin_count(), scan_odim.ray_count()});
            data_odim.read_unpack(scan.data.data(), undetect, nodata);

            vol.sweeps.push_back(std::move(scan));
            break;
        }
    }

    return vol;
}


auto leroi(
    std::filesystem::path const& input_file,
    std::filesystem::path const& out_path,
    io::configuration const& config
) -> void{

    auto proj = map_projection{map_projection::default_context, config["proj4"]};
    auto coords = grid_coordinates{config["size"], config["left_top"], config["cell_delta"], config["units"], config["units"]};
    auto latlons = determine_geodetic_coordinates(proj, coords);
    auto altitudes = init_altitudes(config);
    string moment = config["moment"];

    // Speckle Filter.
    auto min_neighbours = int(config["speckle_min_neighbours"]);
    auto speckle_iters = int(config["speckle_iterations"]);

    std::cout << "Initialisation finished." << std::endl;
    auto vol = read_volume(input_file, config["moment"]);
    std::cout << "Volume read." << std::endl;

    // Setup output product
    auto out_file = io::nc::file{out_path, io_mode::create};
    auto& dim_a = out_file.create_dimension("alt", altitudes.size());
    auto& dim_y = io::cf::create_spatial_dimension(out_file, "y", "projection_y_coordinate", coords.row_units(), coords.row_edges());
    auto& dim_x = io::cf::create_spatial_dimension(out_file, "x", "projection_x_coordinate", coords.col_units(), coords.col_edges());
    auto& var_a = out_file.create_variable("alt", io::nc::data_type::f32, {&dim_a});
    var_a.att_set("standard_name", "altitude");
    var_a.att_set("units", "m");
    var_a.write(altitudes);

    std::cout << "Netcdf initialised." << std::endl;

    auto& var_b = out_file.create_variable(moment, io::nc::data_type::f32, {&dim_a, &dim_y, &dim_x}, {1, dim_y.size(), dim_x.size()});

    for(size_t i = 0; i < altitudes.size(); i++){
        auto cappi = generate_cappi(vol, latlons, config["max_alt_dist"], config["idw_pwr"], (float) altitudes[i]);
        var_b.write(cappi, {i});
    }
    var_b.att_set("units", "dBZ");

    std::cout << "Finished." << std::endl;
}


int main(int argc, char* argv[]){
    try{
        // process command line
        while (true){
            int option_index = 0;
            int c = getopt_long(argc, argv, short_options, long_options, &option_index);
            if (c == -1)
                break;
            switch (c){
                case 'h':
                    std::cout << usage_string;
                    return EXIT_SUCCESS;
                case 'g':
                    std::cout << example_config;
                    return EXIT_SUCCESS;
                case '?':
                    std::cerr << try_again;
                    return EXIT_FAILURE;
                }
        }

        if (argc - optind != 3){
            std::cerr << "missing required parameter\n" << try_again;
            return EXIT_FAILURE;
        }

        std::cout << "Reading file: " << argv[optind + 0] << "\nOutput file: " << argv[optind + 1] << std::endl;

        leroi(
            argv[optind + 0],
            argv[optind + 1],
            io::configuration{std::ifstream{argv[optind + 2]}}
        );
    }
    catch (std::exception& err)
    {
        trace::error("fatal exception: {}", format_exception(err));
        return EXIT_FAILURE;
    }
    catch (...)
    {
        trace::error("fatal exception: (unknown exception)");
        return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
