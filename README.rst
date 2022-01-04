A SunPy-affiliated package for solar instrument-specific tools.
---------------------------------------------------------------

.. image:: http://img.shields.io/badge/powered%20by-SunPy-orange.svg?style=flat
    :target: http://www.sunpy.org
    :alt: Powered by SunPy Badge

What is sunkit-instruments?
---------------------------

sunkit-instruments is a SunPy-affiliated package for solar instrument-specific tools.
Its purpose is not to be a repository for all tools for all instruments.
Instead it is intended to perform three main roles:

1. Hold instrument tools that are so few they do not warrant their own package;
2. Hold tools for instruments with no instrument team or the instrument team does not currently support solar applications;
3. Act as an incubator for instrument-specific tools that can evolve into a separate instrument package, backed by an instrument team.

For instrument teams, this package can act as a forum to engage with their user base and learn how to best develop their tools within the SunPy/scientific Python ecosystem.
It also lowers the barrier to publishing Python-based instrument tools by providing packaging and release infrastructure and support.
However should instrument teams want to develop at their own pace or provide a large number of tools,
they should consider starting their own package for full control.
We encourage and support instrument teams in choosing this route and hope they will still engage and collaborate with the SunPy and wider community during their development.
We point to the recent development of `aiapy <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy>`__ as a great example of this type of collaboration.

License
-------

This project is Copyright (c) The SunPy Developers and licensed under the terms of the BSD 3-Clause license.
This package is based upon the `Openastronomy packaging guide <https://github.com/OpenAstronomy/packaging-guide>`_ which is licensed under the BSD 3-clause licence. See the licenses folder for more information.
