"""
====================================
Plotting NuSTAR Spectral Information
====================================

This example shows how to read NuSTAR spectral files while inspecting their contents.

NuSTAR spectral FITS files are recognized by the filename extension, i.e. the sprectrum file has ``.pha``, the effective area has ``.arf``, and the instrument response has ``.rmf``.

"""

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from parfive import Downloader

import astropy.units as u

from sunkit_instruments.nustar.nustar import NustarSpectrum

###############################################################################
# We start with getting the data. This is done by downloading the data
# from a preprocessed NuSTAR observation.
#
# In this case, we will use requests as to keep this example self
# contained but using your browser will also work.
#
# Using the urls:
# http://foxsi.space.umn.edu/data/tmp/jdhsr/test_nustar_data/nu20619003001A06_cl_grade0_sr.pha
# http://foxsi.space.umn.edu/data/tmp/jdhsr/test_nustar_data/nu20619003001A06_cl_grade0_sr.arf
# http://foxsi.space.umn.edu/data/tmp/jdhsr/test_nustar_data/nu20619003001A06_cl_grade0_sr.rmf

urls = [
    "http://foxsi.space.umn.edu/data/tmp/jdhsr/test_nustar_data/nu20619003001A06_cl_grade0_sr.pha",
    "http://foxsi.space.umn.edu/data/tmp/jdhsr/test_nustar_data/nu20619003001A06_cl_grade0_sr.arf",
    "http://foxsi.space.umn.edu/data/tmp/jdhsr/test_nustar_data/nu20619003001A06_cl_grade0_sr.rmf",
]
files = Downloader.simple_download(urls)

###############################################################################
# Loading in spectral files
# -------------------------
# First let's read in the NuSTAR spectral files into the ``NustarSpectrum``
# class.
#
# This class helps maintain cohesion between the information from all
# three files and assigns the correct units for the spectral files.

nustar_spectrum = NustarSpectrum(pha_file=files[0],
                                 arf_file=files[1],
                                 rmf_file=files[2])

###############################################################################
# The information from each file can be extracted simply.

# get PHA/spectrum information
pha_dict = nustar_spectrum.get_pha_info()
# get ARF/effective area information
arf_dict = nustar_spectrum.get_arf_info()
# get RMF/instrument response information
rmf_dict = nustar_spectrum.get_rmf_info()

###############################################################################
# The spectral response matrix (ARF and RMF combined) can also be obtained.

# get the SRM information
srm_dict = nustar_spectrum.get_srm_info()

###############################################################################
# The PHA and ARF information share different axes with the RMF and SRM
# information. This means that if one is changed, then at least one
# other product's information needs to be changed.
#
# Let's plot what all of these look like.

fig = plt.figure(figsize=(11, 9))
gs = gridspec.GridSpec(2, 2)
# define figure structure
pha_axes = fig.add_subplot(gs[0, 0])
arf_axes = fig.add_subplot(gs[0, 1])
rmf_axes = fig.add_subplot(gs[1, 0])
srm_axes = fig.add_subplot(gs[1, 1])
# plot the spectrum/PHA information
rate = pha_dict["spectrum_counts"]/pha_dict["effective_exposure"]
rate_edges_1d = np.append(pha_dict["spectrum_counts_axis"][:, 0], pha_dict["spectrum_counts_axis"][-1, 1]).value
pha_axes.stairs(rate, rate_edges_1d)
pha_axes.set_xlabel(f"Count/Output Energy [{pha_dict["spectrum_counts_axis"].unit:latex}]")
pha_axes.set_ylabel(f"Rates [{rate.unit:latex}]")
pha_axes.set_title("Observed Spectrum")
pha_axes.set_xscale("log")
pha_axes.set_yscale("log")
# plot the effective area/ARF information
arf_edges_1d = np.append(arf_dict["effective_area_axis"][:, 0], arf_dict["effective_area_axis"][-1, 1]).value
arf_axes.stairs(arf_dict["effective_area"], arf_edges_1d)
arf_axes.set_xlabel(f"Photon/Input Energy [{arf_dict["effective_area_axis"].unit:latex}]")
arf_axes.set_ylabel(f"Effective Area [{arf_dict["effective_area"].unit:latex}]")
arf_axes.set_title("Ancillary Response")
arf_axes.set_xscale("log")
arf_axes.set_yscale("log")
# plot the redistribution/RMF information
rmf_output_edges_1d = np.append(rmf_dict["rmf_output_axis"][:, 0], rmf_dict["rmf_output_axis"][-1, 1]).value
rmf_input_edges_1d = np.append(rmf_dict["rmf_input_axis"][:, 0], rmf_dict["rmf_input_axis"][-1, 1]).value
r = rmf_axes.pcolormesh(rmf_output_edges_1d, rmf_input_edges_1d, rmf_dict["rmf"].value, norm=LogNorm())
cbarr = plt.colorbar(r)
cbarr.ax.set_ylabel(f"Detector Response [{rmf_dict["rmf"].unit:latex}]")
rmf_axes.set_xlabel(f"Count/Output Energy [{rmf_dict["rmf_output_axis"].unit:latex}]")
rmf_axes.set_ylabel(f"Photon/Input Energy [{rmf_dict["rmf_input_axis"].unit:latex}]")
rmf_axes.set_title("Redistribution Matrix")
# plot the spectral response/SRM information
srm_output_edges_1d = np.append(srm_dict["srm_output_axis"][:, 0], srm_dict["srm_output_axis"][-1, 1]).value
srm_input_edges_1d = np.append(srm_dict["srm_input_axis"][:, 0], srm_dict["srm_input_axis"][-1, 1]).value
s = srm_axes.pcolormesh(srm_output_edges_1d, srm_input_edges_1d, srm_dict["srm"].value, norm=LogNorm())
cbars = plt.colorbar(s)
cbars.ax.set_ylabel(f"Spectral Response [{srm_dict["srm"].unit:latex}]")
srm_axes.set_xlabel(f"Count/Output Energy [{srm_dict["srm_output_axis"].unit:latex}]")
srm_axes.set_ylabel(f"Photon/Input Energy [{srm_dict["srm_input_axis"].unit:latex}]")
srm_axes.set_title("Spectral Response Matrix")
# general plot stuff
plt.suptitle("NuSTAR Spectral Information")
plt.tight_layout()
plt.show()

###############################################################################
# The spectrum object
# -------------------
#
# The ``NustarSpectrum`` class will also package the spectral information
# up into a spectrum object.
#
# This object will include the data errors and everything a user needs
# to fit the spectrum in something like `sunkit-spex <https://sunkit-spex.readthedocs.io/en/latest/>`__.

spectrum_object = nustar_spectrum.get_spec_obj()
print(spectrum_object)

###############################################################################
# Advanced: Re-binning
# -------------------
#
# **This should be used with extreme caution and only if a user
# understands what they are doing.**
#
# The `NustarSpectrum` class allows the opportunity to re-bin along the
# the different spectral axes.
#
# To re-bin the data and inspect what the result will be, a user can use
# the functions with the ``rebin_<source>_info`` format where
# ``<source>`` refers to "pha", "arf", "rmf", or "srm".
#
# The user will be warned about the dangers of doing this activity
# separately.
#
# To re-bin the data and apply it to the stored values of the class, we
# can define the new bins we want to use and the ``rebin_info`` method.

# choose new input bins to have 5 keV binning
new_input_axis_edges_1d = np.arange(np.min(rmf_dict["rmf_input_axis"].value),
                                    np.max(rmf_dict["rmf_input_axis"].value),
                                    5) << u.keV
new_input_axis_edges = np.hstack((new_input_axis_edges_1d[:-1,None],
                                   new_input_axis_edges_1d[1:,None])) << u.keV
# choose new output bins to have 1 keV binning
new_output_axis_edges_1d = np.arange(np.min(rmf_dict["rmf_output_axis"].value),
                                     np.max(rmf_dict["rmf_output_axis"].value),
                                     1) << u.keV
new_output_axis_edges = np.hstack((new_output_axis_edges_1d[:-1,None],
                                   new_output_axis_edges_1d[1:,None])) << u.keV
nustar_spectrum.rebin_info(new_input_axis_edges=new_input_axis_edges, new_output_axis_edges=new_output_axis_edges)

###############################################################################
# One or both axes can be updated using the above method.
#
# Inspect the data again to show all axes have been updated appropriately.

# get all the different information components
pha_dict = nustar_spectrum.get_pha_info()
arf_dict = nustar_spectrum.get_arf_info()
rmf_dict = nustar_spectrum.get_rmf_info()
srm_dict = nustar_spectrum.get_srm_info()

fig = plt.figure(figsize=(11, 9))
gs = gridspec.GridSpec(2, 2)
# define figure structure
pha_axes = fig.add_subplot(gs[0, 0])
arf_axes = fig.add_subplot(gs[0, 1])
rmf_axes = fig.add_subplot(gs[1, 0])
srm_axes = fig.add_subplot(gs[1, 1])
# plot the spectrum/PHA information
rate = pha_dict["spectrum_counts"]/pha_dict["effective_exposure"]
rate_edges_1d = np.append(pha_dict["spectrum_counts_axis"][:, 0], pha_dict["spectrum_counts_axis"][-1, 1]).value
pha_axes.stairs(rate, rate_edges_1d)
pha_axes.set_xlabel(f"Count/Output Energy [{pha_dict["spectrum_counts_axis"].unit:latex}]")
pha_axes.set_ylabel(f"Rates [{rate.unit:latex}]")
pha_axes.set_title("Re-binned Observed Spectrum")
pha_axes.set_xscale("log")
pha_axes.set_yscale("log")
# plot the effective area/ARF information
arf_edges_1d = np.append(arf_dict["effective_area_axis"][:, 0], arf_dict["effective_area_axis"][-1, 1]).value
arf_axes.stairs(arf_dict["effective_area"], arf_edges_1d)
arf_axes.set_xlabel(f"Photon/Input Energy [{arf_dict["effective_area_axis"].unit:latex}]")
arf_axes.set_ylabel(f"Effective Area [{arf_dict["effective_area"].unit:latex}]")
arf_axes.set_title("Re-binned Ancillary Response")
arf_axes.set_xscale("log")
arf_axes.set_yscale("log")
# plot the redistribution/RMF information
rmf_output_edges_1d = np.append(rmf_dict["rmf_output_axis"][:, 0], rmf_dict["rmf_output_axis"][-1, 1]).value
rmf_input_edges_1d = np.append(rmf_dict["rmf_input_axis"][:, 0], rmf_dict["rmf_input_axis"][-1, 1]).value
r = rmf_axes.pcolormesh(rmf_output_edges_1d, rmf_input_edges_1d, rmf_dict["rmf"].value, norm=LogNorm())
cbarr = plt.colorbar(r)
cbarr.ax.set_ylabel(f"Detector Response [{rmf_dict["rmf"].unit:latex}]")
rmf_axes.set_xlabel(f"Count/Output Energy [{rmf_dict["rmf_output_axis"].unit:latex}]")
rmf_axes.set_ylabel(f"Photon/Input Energy [{rmf_dict["rmf_input_axis"].unit:latex}]")
rmf_axes.set_title("Re-binned Redistribution Matrix")
# plot the spectral response/SRM information
srm_output_edges_1d = np.append(srm_dict["srm_output_axis"][:, 0], srm_dict["srm_output_axis"][-1, 1]).value
srm_input_edges_1d = np.append(srm_dict["srm_input_axis"][:, 0], srm_dict["srm_input_axis"][-1, 1]).value
s = srm_axes.pcolormesh(srm_output_edges_1d, srm_input_edges_1d, srm_dict["srm"].value, norm=LogNorm())
cbars = plt.colorbar(s)
cbars.ax.set_ylabel(f"Spectral Response [{srm_dict["srm"].unit:latex}]")
srm_axes.set_xlabel(f"Count/Output Energy [{srm_dict["srm_output_axis"].unit:latex}]")
srm_axes.set_ylabel(f"Photon/Input Energy [{srm_dict["srm_input_axis"].unit:latex}]")
srm_axes.set_title("Re-binned Spectral Response Matrix")
# general plot stuff
plt.suptitle("Re-binned NuSTAR Spectral Information")
plt.tight_layout()
plt.show()
