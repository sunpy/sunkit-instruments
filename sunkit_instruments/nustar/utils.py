import numpy as np

import astropy.units as u

__all__ = ["regroup_any_array", "rebin_rmf"]

def regroup_any_array(data:np.ndarray|u.Quantity, old_bins:np.ndarray|u.Quantity, new_bins:np.ndarray|u.Quantity, combine_by:str|None=None):
    """Takes any array of data in old_bins space and rebins along data
    array axis==0 to have new_bins.

    Can specify how the bins are combined: [\"sum\", \"mean\", \"quadrature\"].

    Parameters
    ----------
    data, old_bins, new_bins : `~numpy.ndarray`
            Array of the data, current bins for the data, and new bins
            for the data. Shape of the bin arrays should be `(N+1, 2)`
            where `N` is the length of `data`.

    combine_by : string
            Defines how to combine multiple bins along axis 0. E.g., \"sum\"
            adds the data, \"mean\" averages the data, and \"quadrature\"
            sums the data in quadrature. If `None`, then \"sum\" is used.
            Default: None

    Returns
    -------
    The grouped data array.
    """
    # sanatise inputs
    combine_by = "sum" if combine_by is None else combine_by
    data, du = _get_val_and_unit(data)
    old_bins, obu = _get_val_and_unit(old_bins)
    new_bins, nbu  = _get_val_and_unit(new_bins)
    old_bins = _convert_old_value_to_new_unit_values(old_bins, obu, nbu)

    new_binned_data = []
    for nb in new_bins:
        # just loop through new bins and bin data from between new_bin_lower<=old_bin_lowers and old_bin_highers<new_bin_higher
        if combine_by == "sum":
            new_binned_data.append(
                np.sum(data[np.nonzero((nb[0] <= old_bins[:, 0]) & (nb[-1] >= old_bins[:, -1]))], axis=0)
            )
        elif combine_by == "mean":
            new_binned_data.append(
                np.mean(data[np.nonzero((nb[0] <= old_bins[:, 0]) & (nb[-1] >= old_bins[:, -1]))], axis=0)
            )
        elif combine_by == "quadrature":
            new_binned_data.append(
                np.sqrt(np.sum(data[np.nonzero((nb[0] <= old_bins[:, 0]) & (nb[-1] >= old_bins[:, -1]))] ** 2, axis=0))
            )
    return np.array(new_binned_data) if du is None else np.array(new_binned_data) << du

def rebin_rmf(
    matrix:np.ndarray|u.Quantity, old_output_bins:np.ndarray|u.Quantity=None, new_output_bins:np.ndarray|u.Quantity=None, old_input_bins:np.ndarray|u.Quantity=None, new_input_bins:np.ndarray|u.Quantity=None
):
    """Rebins the photon and/or count channels of the redistribution matrix
    if needed.

    This will rebin any 2d array by taking the mean across the input axis
    (rows or photon space) and summing across the output axis (columns or
    count space).

    Parameters
    ----------
    matrix : 2d array
            Redistribution matrix.

    old_output_bins, new_output_bins : 1d arrays
            The old count channel binning and the new binning to for the
            redistribution matrix columns (sum columns).

    old_input_bins, new_input_bins : 1d arrays
            The old photon channel binning and the new binning to for the
            redistribution matrix columns (average rows).

    Returns
    -------
    The rebinned 2d redistribution matrix.
    """

    matrix, mu = _get_val_and_unit(matrix)
    new_output_bins, nobu = _get_val_and_unit(new_output_bins)
    old_output_bins, oobu = _get_val_and_unit(old_output_bins)
    new_input_bins, nibu = _get_val_and_unit(new_input_bins)
    old_input_bins, oibu = _get_val_and_unit(old_input_bins)

    old_output_bins = _convert_old_value_to_new_unit_values(old_output_bins, oobu, nobu)
    old_input_bins = _convert_old_value_to_new_unit_values(old_input_bins, oibu, nibu)

    # across channel bins, we sum. across energy bins, we average
    # appears to be >2x faster to average first then sum if needing to do both
    if (new_input_bins is not None) and (old_input_bins is not None):
        # very slight difference to rbnrmf when binning across photon axis, <2% of entries have a ratio (my way/rbnrmf) >1 (up to 11)
        # all come from where the original rmf has zeros originally so might be down to precision being worked in, can't expect the exact same numbers essentially
        matrix = regroup_any_array(data=matrix, old_bins=old_input_bins, new_bins=new_input_bins, combine_by="mean")
    if (new_output_bins is not None) and (old_output_bins is not None):
        matrix = regroup_any_array(
            data=matrix.T, old_bins=old_output_bins, new_bins=new_output_bins, combine_by="sum"
        ).T  # need to go along columns so .T then .T back

    return matrix if mu is None else matrix << mu

def _get_val_and_unit(value:np.ndarray|u.Quantity|float|int):
    """Return the value of an object and unit if possible.

    Returns value and `None` if no unit.
    """
    return (value.value, value.unit) if isinstance(value, u.Quantity) else (value, None)

def _convert_old_value_to_new_unit_values(old_value:np.ndarray|float|int, old_unit:u.core.PrefixUnit|u.core.CompositeUnit, new_unit:u.core.PrefixUnit|u.core.CompositeUnit):
    """Convert the old value to new units, return old value if no units.

    Returns unitless value.
    """
    if (old_unit is not None) and (new_unit is not None):
        return ((old_value<<old_unit)<<new_unit).value
    return old_value
