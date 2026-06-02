# LEROI

LEROI is a small Python package for gridding weather radar PPI data.
It interpolates each sweep horizontally with a radius of influence, estimates
the PPI height surfaces, linearly interpolates those sweep fields onto a
Cartesian height grid, and can package the result as a Py-ART `Grid`.

The public API exposed by `import leroi` is:

- `leroi.mask_invalid_data`: despeckle and mask radar fields before gridding.
- `leroi.leroi_interp`: interpolate one or more radar fields onto a Cartesian grid.
- `leroi.build_pyart_grid`: wrap gridded field dictionaries in a Py-ART `Grid`.

## Installation

Install from a local checkout:

```bash
pip install -e .
```

Runtime dependencies are `numpy`, `arm_pyart`, `astropy`, and `scipy`.

## Quick Start

```python
import numpy as np
import pyart
import leroi

radar = pyart.io.read("path/to/radar-file.nc")

fields = ["reflectivity_horizontal"]
radar = leroi.mask_invalid_data(
    radar,
    fields[0],
    add_to=fields,
    correlation_length=2000,
    min_field=0,
    min_area=50,
)

grid_shape = (41, 301, 301)
grid_limits = ((0.0, 20000.0), (-150000.0, 150000.0), (-150000.0, 150000.0))
coords = [
    np.linspace(lower, upper, n)
    for (lower, upper), n in zip(grid_limits, grid_shape)
]

grid_fields = leroi.leroi_interp(
    radar,
    coords,
    field_names=fields,
    weight_type="Barnes",
    Rc=None,
    k=100,
)

grid = leroi.build_pyart_grid(radar, grid_fields, grid_shape, grid_limits)
```

See `notebook/example.ipynb` for a runnable example that reads a Py-ART
compatible radar file, masks a reflectivity field, interpolates it, builds a
Py-ART grid, and plots a PPI next to a grid slice.

## Quality Control

`mask_invalid_data` mutates the input `Radar` object and returns it. It builds a
smoothed copy of the source field, identifies contiguous objects in each PPI
sweep with `pyart.correct.find_objects`, and masks objects smaller than
`min_area`.

Important parameters:

- `field`: source radar field used to identify valid contiguous echo.
- `add_to`: fields that should receive the new mask. Defaults to `[field]`.
- `correlation_length`: along-ray smoothing length in metres for the boxcar
  smoothing step.
- `min_field`: minimum smoothed field value treated as valid echo.
- `min_area`: minimum object area in square kilometres.
- `return_smooth`: keep the temporary `<field>_smooth` field on the radar when
  `True`.

## Gridding

`leroi_interp` returns a dictionary in Py-ART field format:

```python
{
    "field_name": {
        "data": np.ma.MaskedArray,
        "...": "metadata copied from radar.fields[field_name]",
    }
}
```

The field arrays have shape `(z, y, x)` according to `coords`.

Important parameters:

- `coords`: three one-dimensional coordinate arrays ordered as `(z, y, x)`, in
  metres relative to the radar origin.
- `field_names`: a field name or list of names. If omitted, all radar fields are
  interpolated.
- `gatefilter`: optional Py-ART `GateFilter`; excluded gates are removed before
  interpolation.
- `weight_type`: horizontal weighting method, either `"Barnes"` or `"Cressman"`.
  Barnes weighting is the default.
- `Rc`: radius of influence in metres. If `None`, LEROI estimates one from the
  maximum azimuthal gate spacing across the sweeps.
- `k`: maximum nearest neighbours queried around each grid point. Increase this
  if LEROI warns that valid points are being left out of the radius of influence.
- `ground_elevation`: set the first sweep height to zero when its fixed angle is
  at or below this elevation angle.
- `advection`: optional constant `(u, v)` advection in metres per second used to
  shift gate locations relative to `t0`.
- `t0`: reference time in the same units as `radar.time["data"]`.

Use `build_pyart_grid` when downstream tools expect a Py-ART `Grid`. Its
`grid_shape` and `grid_limits` arguments should match the coordinates used for
`leroi_interp`.

## Notes

- Vertical 90 degree scans are ignored when building PPI interpolation surfaces.
- LEROI warns when the requested grid extends outside the radar coverage enough
  that PPI heights cannot be filled everywhere.

## Reference

Brook, J. P., A. Protat, C. K. Potvin, J. S. Soderholm, and H. McGowan, 2023: The 
Effects of Spatial Interpolation on a Novel, Dual-Doppler 3D Wind Retrieval Technique. 
J. Atmos. Oceanic Technol., 40, 1325–1347, https://doi.org/10.1175/JTECH-D-23-0004.1.
