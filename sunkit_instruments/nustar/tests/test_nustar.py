from sunkit_instruments.nustar.nustar import (
    NustarSpectrum,
    )
from sunkit_instruments.nustar.utils import (
    regroup_any_array,
    )

def test_NustarSpectrum():
    base_file = "/Users/kris/Documents/umnPostdoc/projects/analysis/nustarJan2020/spectralFitting/specStuff20515018001/event/specFitting/3638/nu20515018001B06_cl_grade0_sr"
    pha_file = f"{base_file}.pha"
    arf_file = f"{base_file}.arf"
    rmf_file = f"{base_file}.rmf"
    # return NustarSpectrum(pha_file,
    #                arf_file=arf_file,
    #                rmf_file=rmf_file)

# a = test_NustarSpectrum()