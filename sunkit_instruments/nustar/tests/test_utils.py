import astropy.units as u
import numpy as np

from sunkit_instruments.nustar.utils import (
    regroup_any_array,
    rebin_rmf,
    _get_val_and_unit,
    _convert_old_value_to_new_unit_values,
    )

def test_regroup_any_array_sum():
    """Tests for the regrouping function while summing."""

    # testing sum
    orig_bins0 = np.array([[0, 1], [1, 2], [2, 3], [3, 4]])
    orig_data0 = np.array([10, 20, 30, 40])
    new_bins0 = np.array([[0, 2], [2, 4]])
    expected0 = np.array([30, 70])
    result0 = regroup_any_array(orig_data0, orig_bins0, new_bins0, combine_by="sum")
    # test default `combine_by`
    result0a = regroup_any_array(orig_data0, orig_bins0, new_bins0)
    # test with units
    result0b = regroup_any_array(orig_data0<<u.ct/u.s, orig_bins0<<u.keV, new_bins0<<u.keV)
    # test with different, but convertable units
    result0c = regroup_any_array(orig_data0<<u.ct/u.s, orig_bins0<<u.keV, (new_bins0*1e3)<<u.eV)

    np.testing.assert_allclose(expected0, result0)
    np.testing.assert_allclose(result0, result0a)
    np.testing.assert_allclose(result0<<u.ct/u.s, result0b)
    np.testing.assert_allclose(result0<<u.ct/u.s, result0c)
    
    new_bins1 = np.array([[0, 2], [2, 3], [3, 4]])
    expected1 = np.array([30, 30, 40])
    result1 = regroup_any_array(orig_data0, orig_bins0, new_bins1, combine_by="sum")

    new_bins2 = np.array([[1, 2], [2, 3], [3, 4]])
    expected2 = np.array([20, 30, 40])
    result2 = regroup_any_array(orig_data0, orig_bins0, new_bins2, combine_by="sum")

    np.testing.assert_allclose(expected1, result1)
    np.testing.assert_allclose(expected2, result2)

def test_regroup_any_array_mean():
    """Tests for the regrouping function while averaging."""

    # testing mean
    orig_bins0 = np.array([[0, 1], [1, 2], [2, 3], [3, 4]])
    orig_data0 = np.array([10, 20, 30, 40])
    new_bins0 = np.array([[0, 2], [2, 4]])
    expected0 = np.array([15, 35])
    result0 = regroup_any_array(orig_data0, orig_bins0, new_bins0, combine_by="mean")

    new_bins1 = np.array([[0, 2], [2, 3], [3, 4]])
    expected1 = np.array([15, 30, 40])
    result1 = regroup_any_array(orig_data0, orig_bins0, new_bins1, combine_by="mean")

    new_bins2 = np.array([[1, 2], [2, 3], [3, 4]])
    expected2 = np.array([20, 30, 40])
    result2 = regroup_any_array(orig_data0, orig_bins0, new_bins2, combine_by="mean")

    np.testing.assert_allclose(expected0, result0)
    np.testing.assert_allclose(expected1, result1)
    np.testing.assert_allclose(expected2, result2)

def test_regroup_any_array_quadrature():
    """Tests for the regrouping function while averaging."""

    # testing quadrature
    orig_bins0 = np.array([[0, 1], [1, 2], [2, 3], [3, 4]])
    orig_data0 = np.array([10, 20, 30, 40])
    new_bins0 = np.array([[0, 2], [2, 4]])
    expected0 = np.array([np.sqrt(10**2+20**2), np.sqrt(30**2+40**2)])
    result0 = regroup_any_array(orig_data0, orig_bins0, new_bins0, combine_by="quadrature")

    new_bins1 = np.array([[0, 2], [2, 3], [3, 4]])
    expected1 = np.array([np.sqrt(10**2+20**2), 30, 40])
    result1 = regroup_any_array(orig_data0, orig_bins0, new_bins1, combine_by="quadrature")

    new_bins2 = np.array([[1, 2], [2, 3], [3, 4]])
    expected2 = np.array([20, 30, 40])
    result2 = regroup_any_array(orig_data0, orig_bins0, new_bins2, combine_by="quadrature")

    np.testing.assert_allclose(expected0, result0)
    np.testing.assert_allclose(expected1, result1)
    np.testing.assert_allclose(expected2, result2)

def test_rebin_rmf():
    """Test the rebin rmf function."""
    # rebin over input
    orig_rmf0 = np.array([[1, 2, 3], 
                          [4, 5, 6], 
                          [7, 8, 9], 
                          [10, 11, 12]])
    orig_input_bins0 = np.array([[0, 1], [1, 2], [2, 3], [3, 4]])
    orig_output_bins0 = np.array([[10, 20], [20, 30], [30, 40]])
    new_input_bins0 = np.array([[0, 2], [2, 4]])
    result0 = rebin_rmf(
        orig_rmf0, 
        old_output_bins=orig_output_bins0, 
        new_output_bins=None, 
        old_input_bins=orig_input_bins0, 
        new_input_bins=new_input_bins0
        )
    expected0 = np.array([[2.5, 3.5, 4.5], 
                          [8.5, 9.5, 10.5]])
    # rebin over output
    new_output_bins0 = np.array([[10, 30], [30, 40]])
    result1 = rebin_rmf(
        orig_rmf0, 
        old_output_bins=orig_output_bins0, 
        new_output_bins=new_output_bins0, 
        old_input_bins=orig_input_bins0, 
        new_input_bins=None
        )
    expected1 = np.array([[3, 3], 
                          [9, 6], 
                          [15, 9], 
                          [21, 12]])
    # rebin over both axes
    result2 = rebin_rmf(
        orig_rmf0, 
        old_output_bins=orig_output_bins0, 
        new_output_bins=new_output_bins0, 
        old_input_bins=orig_input_bins0, 
        new_input_bins=new_input_bins0
        )
    expected2 = np.array([[6, 4.5],  
                          [18, 10.5]])
    # rebin with units
    result3 = rebin_rmf(
        orig_rmf0<<u.ct/u.ph, 
        old_output_bins=orig_output_bins0<<u.keV, 
        new_output_bins=new_output_bins0<<u.keV, 
        old_input_bins=orig_input_bins0<<u.keV, 
        new_input_bins=new_input_bins0<<u.keV
        )
    expected3 = np.array([[6, 4.5],  
                          [18, 10.5]])<<u.ct/u.ph
    # rebin with different units
    result4 = rebin_rmf(
        orig_rmf0<<u.ct/u.ph, 
        old_output_bins=orig_output_bins0<<u.keV, 
        new_output_bins=(new_output_bins0*1000)<<u.eV, 
        old_input_bins=orig_input_bins0<<u.keV, 
        new_input_bins=new_input_bins0<<u.keV
        )
    expected4 = np.array([[6, 4.5],  
                          [18, 10.5]])<<u.ct/u.ph

    np.testing.assert_allclose(expected0, result0)
    np.testing.assert_allclose(expected1, result1)
    np.testing.assert_allclose(expected2, result2)
    np.testing.assert_allclose(expected3, result3)
    np.testing.assert_allclose(expected4, result4)

def test__get_val_and_unit():
    """Test for `_get_val_and_unit`."""
    a0 = 9
    v0, u0 = _get_val_and_unit(a0)
    assert v0==a0
    assert u0 is None

    a1 = 5 << u.keV
    v1, u1 = _get_val_and_unit(a1)
    assert v1==a1.value
    assert u1 is u.keV

    a2 = np.array([3, 8]) << u.keV
    v2, u2 = _get_val_and_unit(a2)
    np.testing.assert_allclose(v2, a2.value)
    assert u2 is u.keV

def test__convert_old_value_to_new_unit_values():
    """Test for `_convert_old_value_to_new_unit_values`."""
    old_value, old_unit, new_unit = 5, u.keV, u.eV
    output = _convert_old_value_to_new_unit_values(old_value, old_unit, new_unit)
    assert output==5000
