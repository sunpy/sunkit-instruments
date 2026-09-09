import numpy as np

def regroup_any_array(data, old_bins, new_bins, combine_by="sum"):
    """Takes any array of data in old_bins space and rebins along data array axis==0 to have new_bins.

    Can specify how the bins are combined.

    Parameters
    ----------
    data, old_bins, new_bins : np.array
            Array of the data, current bins (for data axis==0), and new bins (for data axis==0).
            Need len(data)==len(old_bins).

    combine_by : string
            Defines how to combine multiple bins along axis 0. E.g., "sum" adds the data, "mean" averages
            the data, and "quadrature" sums the data in quadrature.
            Default: "sum"

    Returns
    -------
    The new bins and the corresponding grouped counts.
    """
    new_binned_data = []
    for nb in new_bins:
        # just loop through new bins and bin data from between new_bin_lower<=old_bin_lowers and old_bin_highers<new_bin_higher
        if combine_by == "sum":
            new_binned_data.append(
                np.sum(data[np.where((nb[0] <= old_bins[:, 0]) & (nb[-1] >= old_bins[:, -1]))], axis=0)
            )
        elif combine_by == "mean":
            new_binned_data.append(
                np.mean(data[np.where((nb[0] <= old_bins[:, 0]) & (nb[-1] >= old_bins[:, -1]))], axis=0)
            )
        elif combine_by == "quadrature":
            new_binned_data.append(
                np.sqrt(np.sum(data[np.where((nb[0] <= old_bins[:, 0]) & (nb[-1] >= old_bins[:, -1]))] ** 2, axis=0))
            )
    return np.array(new_binned_data)

def rebin_rmf(
    matrix, old_count_bins=None, new_count_bins=None, old_photon_bins=None, new_photon_bins=None, axis="count"
):
    """Rebins the photon and/or count channels of the redistribution matrix if needed.

    This will rebin any 2d array by taking the mean across photon space (rows) and summing
    across count space (columns).

    If no effective area information from the instrument then this is passed straight
    to `_rebin_srm`, if there is then the `_rebin_srm` should be overwritten.

    Parameters
    ----------
    matrix : 2d array
            Redistribution matrix.

    old_count_bins, new_count_bins : 1d arrays
            The old count channel binning and the new binning to be for the redistribution matrix columns (sum columns).

    old_photon_bins, new_photon_bins : 1d arrays
            The old photon channel binning and the new binning to be for the redistribution matrix columns (average rows).

    axis : string
            Define what \'axis\' the binning should be applied to. E.g., \'photon\', \'count\', or \'photon_and_count\'.
            Default: \'count\'

    Returns
    -------
    The rebinned 2d redistribution matrix.
    """
    # across channel bins, we sum. across energy bins, we average
    # appears to be >2x faster to average first then sum if needing to do both
    if (axis == "photon") or (axis == "photon_and_count"):
        # very slight difference to rbnrmf when binning across photon axis, <2% of entries have a ratio (my way/rbnrmf) >1 (up to 11)
        # all come from where the original rmf has zeros originally so might be down to precision being worked in, can't expect the exact same numbers essentially
        matrix = regroup_any_array(data=matrix, old_bins=old_photon_bins, new_bins=new_photon_bins, combine_by="mean")
    if (axis == "count") or (axis == "photon_and_count"):
        matrix = regroup_any_array(
            data=matrix.T, old_bins=old_count_bins, new_bins=new_count_bins, combine_by="sum"
        ).T  # need to go along columns so .T then .T back

    return matrix
