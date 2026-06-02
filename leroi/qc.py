import pyart
import numpy as np
from astropy.convolution import convolve
import warnings

def smooth_ppi(radar, field, sweep, c_len):
    """Return a smoothed PPI sweep for a radar field.

    Smoothing is performed along each ray with a one-dimensional boxcar kernel.
    `astropy.convolution.convolve` fills short NaN gaps during convolution, and
    `gate_range_mask` masks edge gates where the smoothing window is not fully
    supported by valid data.

    Parameters
    ----------
    radar : pyart.core.Radar
        Radar object containing the source field.
    field : str
        Radar field name.
    sweep : int
        Sweep number to smooth.
    c_len : float
        Smoothing length in metres.

    Returns
    -------
    np.ma.MaskedArray
        Smoothed two-dimensional PPI data for the requested sweep.
    """
    window = int(np.ceil(c_len / np.mean(np.diff(radar.range["data"]))) // 2 * 2 + 1)
    data = radar.get_field(sweep, field).filled(np.nan)
    kernel = np.ones((1, window)) / float(window)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        smooth = convolve(data, kernel, boundary="extend")
    mask = gate_range_mask(smooth, window)
    return np.ma.masked_array(smooth, np.logical_or(np.isnan(smooth), mask))


def gate_range_mask(data, window):
    """Return a mask for unsupported gates at the edges of each ray.

    Parameters
    ----------
    data : array-like
        Two-dimensional PPI data with shape `(rays, gates)`.
    window : int
        Smoothing window length in gates.

    Returns
    -------
    np.ndarray
        Boolean mask with `True` where the smoothed data should be masked.
    """
    window -= int(window / 2)
    end = data.shape[1]
    mask = ~np.isnan(data)
    starts = np.argmax(mask, axis=1)
    starts[starts > window] += window - 1
    ends = np.argmax(mask[:, ::-1], axis=1)
    ends = end - ends
    ends[ends < (end - window)] += 1 - window
    mask = np.ones(data.shape)
    for j, i in enumerate(zip(starts, ends)):
        mask[j, i[0] : i[1]] = 0
    return mask.astype("bool")


def _clear_small_echoes_ppi(label_image, areas, min_area):
    """Remove labelled PPI objects with total area below `min_area`.

    Parameters
    ----------
    label_image : np.ndarray
        Integer object labels, as returned by `pyart.correct.find_objects`.
    areas : np.ndarray
        Gate area estimates in square kilometres.
    min_area : float
        Minimum object area in square kilometres.

    Returns
    -------
    np.ndarray
        Label image with small objects set to zero.
    """
    small_echoes = []
    for i in range(1, np.amax(label_image) + 1):
        area = np.sum(areas[label_image == i])
        if area < min_area:
            small_echoes.append(i)
    small = np.array(small_echoes)
    for obj in small:
        label_image[label_image == obj] = 0
        label_image[label_image > obj] -= 1
        small[small > obj] -= 1
    return label_image


def mask_invalid_data(radar, field, add_to=None, correlation_length=2000, 
                      min_field=0, min_area=10, return_smooth=False):
    """Mask radar data outside contiguous PPI echo objects.

    The function smooths `field`, identifies contiguous objects in each non-
    vertical PPI sweep with `pyart.correct.find_objects`, removes objects whose
    area is below `min_area`, and applies the resulting mask to `add_to` fields.
    The input radar is modified in place and returned.

    Parameters
    ----------
    radar : pyart.core.Radar
        Radar object to modify.
    field : str
        Field used to identify contiguous valid echo.
    add_to : list[str] or None, optional
        Field names that should receive the new mask. Defaults to `[field]`.
    correlation_length : float, optional
        Along-ray smoothing length in metres.
    min_field : float, optional
        Minimum smoothed field value treated as valid echo.
    min_area : float, optional
        Minimum contiguous object area in square kilometres.
    return_smooth : bool, optional
        Keep the temporary `<field>_smooth` field when `True`.

    Returns
    -------
    pyart.core.Radar
        The same radar object with updated field masks.
    """
    dR = np.mean(np.diff(radar.range["data"])) / 1e3
    elevations = radar.fixed_angle['data']
    smooth_data, mask = np.ma.zeros((2, radar.nrays, radar.ngates))
    smooth_fn = field + "_smooth"
    nsweeps = radar.nsweeps
    for sweep in range(nsweeps):
        smooth_data[radar.get_slice(sweep)] = smooth_ppi(radar, field, sweep, correlation_length)
    radar.add_field_like(field, smooth_fn, smooth_data, replace_existing=True)
    for sweep in range(nsweeps):
        if elevations[sweep] == 90.0:
            continue # dont mask vertical scans
        dA = np.radians(np.mean(np.diff(np.sort(radar.azimuth["data"][radar.get_slice(sweep)]))))
        ranges = np.repeat(radar.range["data"][np.newaxis, :], radar.rays_per_sweep["data"][sweep], axis=0)
        corr = np.cos(np.radians(radar.fixed_angle["data"][sweep])) ** 2
        areas = ranges * dA * dR * corr / 1e3
        obj_dict = pyart.correct.find_objects(radar, smooth_fn, min_field, sweeps=sweep)
        label = _clear_small_echoes_ppi(obj_dict["data"].filled(0), areas, min_area)
        mask[radar.get_slice(sweep)] = label == 0
    if add_to is None:
        add_to = [
            field,
        ]
    for f in add_to:
        radar.fields[f]["data"] = np.ma.masked_array(radar.fields[f]["data"], mask)
    if not return_smooth:
        radar.fields.pop(smooth_fn)
    return radar
