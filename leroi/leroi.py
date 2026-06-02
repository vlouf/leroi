"""Interpolation utilities for gridding Py-ART radar PPI data."""

import time
import warnings
import multiprocessing as mp

import pyart
import numpy as np
from scipy.spatial import cKDTree
from astropy.convolution import convolve
from pyart.config import get_metadata

def fill_heights(a):
    """Fill missing PPI heights above or below valid data in each grid column.

    Missing heights below the lowest valid sweep are filled with that lowest
    valid height, and missing heights above the highest valid sweep are filled
    with the highest valid height. Remaining missing values indicate columns
    that are outside radar coverage or have gaps inside the sweep stack, and a
    warning is issued.

    Parameters
    ----------
    a : np.ma.MaskedArray
        PPI height array with shape `(sweep, y, x)`.

    Returns
    -------
    np.ma.MaskedArray
        Height array with fillable edge gaps replaced.
    """
    valid = ~a.mask
    a = a.filled(np.nan)
    bidx = valid.argmax(axis=0)
    tidx = valid.shape[0]-np.flip(valid, axis = 0).argmax(axis=0)
    J,I = np.where(bidx!=0)
    for j,i in zip(J,I):
        a[:bidx[j,i],j,i] = a[bidx[j,i],j,i]
    J,I = np.where(tidx!=3)
    for j,i in zip(J,I):
        a[tidx[j,i]:,j,i] = a[(tidx[j,i]-1),j,i]
    out = np.ma.masked_invalid(a)
    if (out.mask).sum() >0:
        warnings.warn("""There are still missing heights despite filling.
        This means the grid is OUTSIDE the radar range!""")
    return out

def get_data_mask(radar, fields, gatefilter=None):
    """Create a combined invalid-data mask for interpolation.

    A gate is valid only when it is unmasked in every requested field and, when
    supplied, is not excluded by the Py-ART gate filter.

    Parameters
    ----------
    radar : pyart.core.Radar
        Radar object containing the fields.
    fields : list[str]
        Field names whose masks should be combined.
    gatefilter : pyart.correct.GateFilter or None, optional
        Optional Py-ART gate filter.

    Returns
    -------
    np.ndarray
        Boolean mask with the same `(ray, gate)` shape as the radar fields.
        `True` marks invalid gates.
    """

    mask = np.ones((radar.nrays, radar.ngates)).astype("bool")

    # combine data masks
    for field in fields:
        mask *= ~radar.fields[field]["data"].mask

    # combine gatefilter
    if gatefilter is not None:
        mask *= ~gatefilter.gate_excluded
    return ~mask


def get_leroy_roi(radar, coords, frac=0.6):
    """Estimate the horizontal radius of influence from sweep azimuth spacing.

    The estimate follows the spacing guidance used by Dahl et al. (2019): LEROI
    finds the largest azimuthal spacing among PPI sweeps and scales the maximum
    grid range by `frac`.

    Parameters
    ----------
    radar : pyart.core.Radar
        Radar object.
    coords : sequence[np.ndarray]
        One-dimensional grid coordinates ordered as `(z, y, x)`.
    frac : float, optional
        Fraction of the largest azimuthal spacing to use.

    Returns
    -------
    float
        Radius of influence in metres.
    """
    roi = 0
    rmax = np.sqrt(max(abs(coords[0])) ** 2 + max(abs(coords[1])) ** 2 + max(abs(coords[2])) ** 2)
    sort_idx = np.argsort(radar.fixed_angle["data"])
    for i in sort_idx:
        az = np.amax(np.radians(np.amax(np.diff(np.sort(radar.azimuth["data"][radar.get_slice(i)])))))
        r = frac * az * rmax
        if r > roi:
            roi = r
    return roi


def _calculate_ppi_heights(radar, coords, weight_type, Rc, multiprocessing, ground_elevation):
    """Calculate PPI sweep heights on the horizontal grid.

    Each sweep's gate heights are projected onto the requested `(y, x)` grid
    with the same Barnes or Cressman radius-of-influence weighting used for
    fields.

    Parameters
    ----------
    radar : pyart.core.Radar
        Radar object.
    coords : sequence[np.ndarray]
        One-dimensional grid coordinates ordered as `(z, y, x)`.
    weight_type : {"Barnes", "Cressman"}
        Horizontal weighting method.
    Rc : float
        Radius of influence in metres.
    multiprocessing : bool
        Use all available `scipy.spatial.cKDTree.query` workers when `True`.
    ground_elevation : float
        Set the lowest sweep height to zero when its fixed angle is at or below
        this value.

    Returns
    -------
    np.ma.MaskedArray
        Sweep height array with shape `(sweep, y, x)`.
    """
    slices = []
    elevations = radar.fixed_angle["data"]

    # multiprocessing, 1 for none, -1 for all available workers
    workers = (-1) ** int(multiprocessing)

    # sort sweep index to process from lowest sweep and ascend
    sort_idx = list(np.argsort(elevations))
    if 90.0 in elevations:
        sort_idx.remove(np.argwhere(elevations == 90))
    Y, X = np.meshgrid(coords[1], coords[2], indexing="ij")
    for j,i in enumerate(sort_idx):
        x, y, z = radar.get_gate_x_y_z(i)
        data = z.ravel()
        tree = cKDTree(np.c_[y.ravel(), x.ravel()])
        d, idx = tree.query(np.c_[Y.ravel(), X.ravel()], k=10, distance_upper_bound=Rc, workers=workers)
        idx[idx == len(data)] = 0

        # do all of the weighting stuff based on kdtree distance
        d[np.isinf(d)] = Rc + 1e5
        d2, r2 = d ** 2, Rc ** 2
        if weight_type == 'Barnes':
            w = np.exp(-(d2) / (r2/4))
            w[w<np.exp(-4)] = 0
        elif weight_type == 'Cressman':
            w = (r2 - d2) / (r2 + d2)
            w[w < 0] = 0
        else:
            raise NotImplementedError("'weight_type' must be either 'Barnes' or 'Cressman'")
        
        sw = np.sum(w, axis=1)
        valid = sw != 0
        # put valid data into a resultant array and reshape to model grid
        slce = np.zeros(sw.shape)
        # set lowest scan heights to zeros if requested
        if (j == 0) and (elevations[i] <= ground_elevation):
            pass
        # dont worry if there's no valid data, masked anyway
        elif len(data) == 0:
            pass
        # perform weighted interp on valid heights
        else:
            slce[valid] = np.sum(data[idx] * w, axis=1)[valid] / sw[valid]

        slce = np.ma.masked_array(slce, mask=~valid)
        slices.append(slce.reshape((len(coords[1]), len(coords[2]))))

    return np.ma.stack(slices)


def smooth_grid(grid, coords, verbose, kernel=None, corr_lens=None, 
                filter_its=0, preserve_nan=False, boundary = "extend"):
    """Smooth a gridded field as a postprocessing step.

    Provide either a custom Astropy convolution `kernel` or correlation lengths
    from which a three-dimensional boxcar kernel can be derived.

    Parameters
    ----------
    grid : array-like
        Three-dimensional gridded field with shape `(z, y, x)`.
    coords : sequence[np.ndarray]
        One-dimensional grid coordinates ordered as `(z, y, x)`.
    verbose : bool
        Print progress messages when `True`.
    kernel : astropy.convolution.Kernel or None, optional
        Kernel passed to `astropy.convolution.convolve`.
    corr_lens : tuple[float, float] or None, optional
        Vertical and horizontal correlation lengths in metres. Required when
        `kernel` is not supplied and `filter_its` is greater than zero.
    filter_its : int, optional
        Number of smoothing iterations.
    preserve_nan : bool, optional
        Passed through to `astropy.convolution.convolve`.
    boundary : str, optional
        Boundary mode passed through to `astropy.convolution.convolve`.

    Returns
    -------
    np.ndarray
        Smoothed grid.
    """

    dz = np.mean(np.diff(np.sort(coords[0])))
    dy = np.mean(np.diff(np.sort(coords[1])))
    dx = np.mean(np.diff(np.sort(coords[2])))
    dh = np.mean((dy, dx))

    if verbose:
        print("Filtering...")

    if kernel is None:
        if corr_lens == None:
            raise NotImplementedError(
                """You must either input a convolution kernel
                ('kernel') or some correlation lengths ('corr_len')."""
            )
        v_window = int(np.ceil(corr_lens[0] / dz) // 2 * 2 + 1)
        h_window = int(np.ceil(corr_lens[1] / dh) // 2 * 2 + 1)
        kernel = np.ones((v_window, h_window, h_window)) / float(v_window * h_window * h_window)

    smooth = grid.copy()
    for i in range(filter_its):
        smooth = convolve(smooth, kernel, boundary=boundary, preserve_nan=preserve_nan)

    return smooth.copy()


def interp_along_axis(y, x, newx, axis, inverse=False, method="linear"):
    """Interpolate values along one axis of an irregular grid.

    This vectorized routine is adapted from Peter Kalverla's March 2018
    StackOverflow answer:
    https://stackoverflow.com/questions/28934767/best-way-to-interpolate-a-numpy-ndarray-along-an-axis

    Parameters
    ----------
    y : np.ndarray
        Values to interpolate.
    x : np.ndarray
        Existing coordinates for `y`.
    newx : np.ndarray
        Target coordinates.
    axis : int
        Axis along which to interpolate.
    inverse : bool, optional
        Reverse the interpolation axis before and after interpolation. This is
        useful for coordinates that are naturally descending.
    method : {"linear", "cubic"}, optional
        Interpolation method.

    Returns
    -------
    np.ndarray
        Interpolated values with `axis` resized to match `newx`.

    Notes
    -----
    The coordinate arrays must be monotonic along `axis`; target values outside
    the source coordinate range are filled with NaN.
    """
    # View of x and y with axis as first dimension
    if inverse:
        _x = np.moveaxis(x, axis, 0)[::-1, ...]
        _y = np.moveaxis(y, axis, 0)[::-1, ...]
        _newx = np.moveaxis(newx, axis, 0)[::-1, ...]
    else:
        _y = np.moveaxis(y, axis, 0)
        _x = np.moveaxis(x, axis, 0)
        _newx = np.moveaxis(newx, axis, 0)

    # Sanity checks
    if np.any(_newx[0] < _x[0]) or np.any(_newx[-1] > _x[-1]):
        # raise ValueError('This function cannot extrapolate')
        warnings.warn("Some values are outside the interpolation range. " "These will be filled with NaN")
    if np.any(np.diff(_x, axis=0) < 0):
        raise ValueError("x should increase monotonically")
    if np.any(np.diff(_newx, axis=0) < 0):
        raise ValueError("newx should increase monotonically")
    # Cubic interpolation needs the gradient of y in addition to its values
    if method == "cubic":
        # For now, simply use a numpy function to get the derivatives
        # This produces the largest memory overhead of the function and
        # could alternatively be done in passing.
        ydx = np.gradient(_y, axis=0, edge_order=2)
    # This will later be concatenated with a dynamic '0th' index
    ind = [i for i in np.indices(_y.shape[1:])]
    # Allocate the output array
    original_dims = _y.shape
    newdims = list(original_dims)
    newdims[0] = len(_newx)
    newy = np.zeros(newdims)
    # set initial bounds
    i_lower = np.zeros(_x.shape[1:], dtype=int)
    i_upper = np.ones(_x.shape[1:], dtype=int)
    x_lower = _x[0, ...]
    x_upper = _x[1, ...]
    for i, xi in enumerate(_newx):
        # Start at the 'bottom' of the array and work upwards
        # This only works if x and newx increase monotonically
        # Update bounds where necessary and possible
        needs_update = (xi > x_upper) & (i_upper + 1 < len(_x))
        # print x_upper.max(), np.any(needs_update)
        while np.any(needs_update):
            i_lower = np.where(needs_update, i_lower + 1, i_lower)
            i_upper = i_lower + 1
            x_lower = _x[tuple([i_lower] + ind)]
            x_upper = _x[tuple([i_upper] + ind)]
            # Check again
            needs_update = (xi > x_upper) & (i_upper + 1 < len(_x))
        # Express the position of xi relative to its neighbours
        xj = (xi - x_lower) / (x_upper - x_lower)
        # Determine where there is a valid interpolation range
        within_bounds = (_x[0, ...] < xi) & (xi < _x[-1, ...])
        if method == "linear":
            f0, f1 = _y[tuple([i_lower] + ind)], _y[tuple([i_upper] + ind)]
            a = f1 - f0
            b = f0
            newy[i, ...] = np.where(within_bounds, a * xj + b, np.nan)
        elif method == "cubic":
            f0, f1 = _y[tuple([i_lower] + ind)], _y[tuple([i_upper] + ind)]
            df0, df1 = ydx[tuple([i_lower] + ind)], ydx[tuple([i_upper] + ind)]
            a = 2 * f0 - 2 * f1 + df0 + df1
            b = -3 * f0 + 3 * f1 - 2 * df0 - df1
            c = df0
            d = f0
            newy[i, ...] = np.where(within_bounds, a * xj ** 3 + b * xj ** 2 + c * xj + d, np.nan)
        else:
            raise ValueError("invalid interpolation method" "(choose 'linear' or 'cubic')")
    if inverse:
        newy = newy[::-1, ...]
    return np.moveaxis(newy, 0, axis)


def _setup_interpolate(radar, coords, dmask, weight_type, Rc, multiprocessing, k, verbose, advection, t0):
    """Precompute sweep interpolation weights and gate indices.

    This builds one horizontal `cKDTree` lookup per PPI sweep and stores only the
    model grid points that have at least one gate inside the radius of influence.
    """
    # setup stuff
    nsweeps = radar.nsweeps
    weights, idxs = [], []
    elevations = radar.fixed_angle["data"]
    Y, X = np.meshgrid(coords[1], coords[2], indexing="ij")
    trim = 0
    model_idxs = -np.ones((nsweeps, len(coords[1]) * len(coords[2])))
    sws, model_lens = [], []

    # multiprocessing, 1 for none, -1 for all available workers
    workers = (-1) ** int(multiprocessing)

    # sort sweep index to process from lowest sweep and ascend
    sort_idx = list(np.argsort(elevations))
    # ignore birdbaths as they aren't ppis
    if 90.0 in elevations:
        sort_idx.remove(np.argwhere(elevations == 90))

    # loop through ppis
    for j, i in enumerate(sort_idx):
        x, y, z = radar.get_gate_x_y_z(i)
        x -= advection[0] * (radar.time['data'] - t0)[radar.get_slice(i),np.newaxis]
        y -= advection[1] * (radar.time['data'] - t0)[radar.get_slice(i),np.newaxis]
        mask = ~dmask[radar.get_slice(i)].flatten()

        # dont bother with ckdtree if there's no data
        if mask.sum() == 0:
            weights.append(np.empty((1, 1)).astype(int))
            idxs.append(np.empty((1, 1)).astype(int))
            sws.append(np.zeros(len(coords[1]) * len(coords[2])))
            model_lens.append(0)
            continue

        # define a lookup tree for the horizontal coordinates
        valid_radar_points = np.c_[y.ravel()[mask], x.ravel()[mask]]
        ndata = valid_radar_points.shape[0]
        tree = cKDTree(valid_radar_points)
        d, idx = tree.query(np.c_[Y.ravel(), X.ravel()], k=k, distance_upper_bound=Rc, workers=workers)

        # check if any kth weight is valid, and trim if possible
        valid = ~(idx == ndata)
        if np.any(valid):
            kidx = max(np.where(valid == 1)[1]) + 1
        else:
            # no valid points from KDTree (valid data coming are outside all ROI regions (outside grid))
            weights.append(np.empty((1, 1)).astype(int))
            idxs.append(np.empty((1, 1)).astype(int))
            sws.append(np.zeros(len(coords[1]) * len(coords[2])))
            model_lens.append(0)
            continue

        if valid[:, -1].sum() > 0:
            warnings.warn("\n Some points are being left out of radius of influence, make 'k' bigger!")

        # set invalid indices to 0 to avoid errors, they are masked out by the weights anyway
        idx[idx == ndata] = 0

        # do all of the weighting stuff based on kdtree distance
        d[np.isinf(d)] = Rc + 1e5
        d2, r2 = d ** 2, Rc ** 2
        if weight_type == 'Barnes':
            w = np.exp(-(d2) / (r2/4))
            w[w<np.exp(-4)] = 0
        elif weight_type == 'Cressman':
            w = (r2 - d2) / (r2 + d2)
            w[w < 0] = 0
        else:
            raise NotImplementedError("'weight_type' must be either 'Barnes' or 'Cressman'")
            
        sw = np.sum(w, axis=1)
        model_idx = np.where(sw != 0)[0]
        model_idxs[j][: len(model_idx)] = model_idx
        model_lens.append(len(model_idx))
        sws.append(sw)
        weights.append(w[model_idx, :kidx])
        idxs.append(idx[model_idx, :kidx])

    # stack weights
    return weights, idxs, model_idxs[:, : max(model_lens)].astype(int), np.array(sws), model_lens


def leroi_interp(
    radar,
    coords,
    field_names=None,
    gatefilter=None,
    weight_type = 'Barnes',
    Rc=None,
    k=100,
    verbose=True,
    smooth_kw = {'filter_its':0},
    multiprocessing=True,
    ground_elevation=-999,
    advection = (0,0),
    t0 = 0
):
    """Interpolate one or more radar fields onto a Cartesian grid.

    LEROI first interpolates each PPI sweep horizontally with Barnes or Cressman
    weighting inside a radius of influence, then linearly interpolates between
    sweep-height surfaces onto the requested vertical coordinate. The method is
    based on Dahl et al. (2019).

    Parameters
    ----------
    radar : pyart.core.Radar
        Radar object containing the source fields.
    coords : sequence[np.ndarray]
        Three one-dimensional grid coordinate arrays ordered as `(z, y, x)`, in
        metres relative to the radar origin.
    field_names : str, list[str], or None, optional
        Field name or names to interpolate. Defaults to all radar fields.
    gatefilter : pyart.correct.GateFilter or None, optional
        Optional Py-ART gate filter; excluded gates are ignored.
    weight_type : {"Barnes", "Cressman"}, optional
        Horizontal radius-of-influence weighting method.
    Rc : float or None, optional
        Radius of influence in metres. When `None`, LEROI estimates a value from
        sweep azimuth spacing.
    k : int, optional
        Maximum number of nearest gates queried around each grid point.
    verbose : bool, optional
        Print progress messages when `True`.
    smooth_kw : dict, optional
        Keyword arguments for `smooth_grid`. The default disables smoothing.
    multiprocessing : bool, optional
        Use all available `scipy.spatial.cKDTree.query` workers when `True`.
    ground_elevation : float, optional
        Set the first sweep height to zero when its fixed angle is at or below
        this value.
    advection : tuple[float, float], optional
        Constant `(u, v)` advection in metres per second applied to gate
        locations relative to `t0`.
    t0 : float, optional
        Reference time in the same units as `radar.time["data"]`.

    Returns
    -------
    dict
        Py-ART-style field dictionary whose arrays have shape `(z, y, x)`.
    """
    ti = time.time()

    if field_names is None:
        # No field defined. Processing all fields in radar.
        field_names = [*radar.fields.keys()]
    if type(field_names) != list:
        field_names = [field_names]

    fields = {}
    dims = [len(coord) for coord in coords]

    if Rc is None:
        Rc = get_leroy_roi(radar, coords)
        if verbose:
            print("Radius of influence set to {} m.".format(Rc))

    dmask = get_data_mask(radar, field_names, gatefilter)
    weights, idxs, model_idxs, sw, model_lens = _setup_interpolate(
        radar, coords, dmask, weight_type, Rc, multiprocessing, k, verbose, advection, t0
    )
    Z, Y, X = np.meshgrid(coords[0], coords[1], coords[2], indexing="ij")
    ppi_height = _calculate_ppi_heights(radar, coords, weight_type, Rc, multiprocessing, ground_elevation)
    # fill in the nan heights. If there is no height data, there will be no radar data by definition.
    ppi_height = fill_heights(ppi_height)
    
    elevations = radar.fixed_angle["data"]

    # sort sweep index to process from lowest sweep and ascend
    sort_idx = list(np.argsort(elevations))
    if 90.0 in elevations:
        sort_idx.remove(np.argwhere(elevations == 90))

    for field in field_names:
        ppis = np.zeros((radar.nsweeps, dims[1] * dims[2]))
        mask = np.ones((radar.nsweeps, dims[1] * dims[2]))
        for i, j in enumerate(sort_idx):
            slc = radar.get_slice(j)
            data = radar.fields[field]["data"].filled(0)[slc][~dmask[slc]]
            if len(data) > 0:
                ppis[i, model_idxs[i, : model_lens[i]]] = (
                    np.sum(data[idxs[i]] * weights[i], axis=1) / sw[i, model_idxs[i, : model_lens[i]]]
                )
                mask[i, model_idxs[i, : model_lens[i]]] = 0
        out = np.ma.masked_array(
            ppis.reshape((radar.nsweeps, dims[1], dims[2])), mask.reshape((radar.nsweeps, dims[1], dims[2]))
        )
        grid = interp_along_axis(out.filled(np.nan), ppi_height, Z, axis=0, method="linear")

        if smooth_kw['filter_its'] > 0:
            grid = smooth_grid(grid, coords, verbose, **smooth_kw)

        # add to output dictionary
        fields[field] = {"data": np.ma.masked_array(grid, mask=np.isnan(grid))}
        # copy the metadata from the radar to the grid
        for key in radar.fields[field].keys():
            if key == "data":
                continue
            fields[field][key] = radar.fields[field][key]

    if verbose:
        print("Took: ", time.time() - ti)

    return fields


def build_pyart_grid(radar, fields, gs, gb, origin = None):
    """Build a Py-ART `Grid` from gridded LEROI fields.

    Parameters
    ----------
    radar : pyart.core.Radar
        Source radar object used for metadata, time, and radar location.
    fields : dict
        Py-ART-style field dictionary, typically returned by `leroi_interp`.
    gs : tuple[int, int, int]
        Grid shape ordered as `(z, y, x)`.
    gb : tuple[tuple[float, float], tuple[float, float], tuple[float, float]]
        Grid bounds ordered as `((z_min, z_max), (y_min, y_max), (x_min, x_max))`.
    origin : tuple[float, float, float] or None, optional
        Grid origin as `(altitude, latitude, longitude)`. Defaults to the radar
        location.

    Returns
    -------
    pyart.core.Grid
        Grid containing the supplied fields.
    """
    # time dictionaries
    time = get_metadata("grid_time")
    time["data"] = np.array([radar.time["data"][0]])
    time["units"] = radar.time["units"]

    # metadata dictionary
    metadata = dict(radar.metadata)

    # grid origin location dictionaries
    origin_latitude = get_metadata("origin_latitude")
    origin_longitude = get_metadata("origin_longitude")
    origin_altitude = get_metadata("origin_altitude")
    
    if origin is None:
        origin_latitude["data"] = radar.latitude["data"]
        origin_longitude["data"] = radar.longitude["data"]
        origin_altitude["data"] = radar.altitude["data"]
    else:
        origin_altitude["data"]  = [origin[0]]
        origin_latitude["data"]  = [origin[1]]
        origin_longitude["data"] = [origin[2]]

    # grid coordinate dictionaries
    nz, ny, nx = gs
    (z0, z1), (y0, y1), (x0, x1) = gb
    x = get_metadata("x")
    x["data"] = np.linspace(x0, x1, nx)
    y = get_metadata("y")
    y["data"] = np.linspace(y0, y1, ny)
    z = get_metadata("z")
    z["data"] = np.linspace(z0, z1, nz)

    # create radar_ dictionaries
    radar_latitude = get_metadata("radar_latitude")
    radar_latitude["data"] = radar.latitude["data"]
    radar_longitude = get_metadata("radar_longitude")
    radar_longitude["data"] = radar.longitude["data"]
    radar_altitude = get_metadata("radar_altitude")
    radar_altitude["data"] = radar.altitude["data"]
    radar_time = time
    radar_name = get_metadata("radar_name")
    name_key = "instrument_name"
    radar_name["data"] = radar.metadata[name_key]
    if name_key in radar.metadata:
        radar_name["data"] = np.array([radar.metadata[name_key]])
    else:
        radar_name["data"] = np.array([""])

    return pyart.core.Grid(
        time,
        fields,
        metadata,
        origin_latitude,
        origin_longitude,
        origin_altitude,
        x,
        y,
        z,
        radar_latitude=radar_latitude,
        radar_longitude=radar_longitude,
        radar_altitude=radar_altitude,
        radar_name=radar_name,
        radar_time=radar_time,
        projection=None,
    )
