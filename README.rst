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
We point to the recent development of `aiapy <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy>`__ as a shining example of this type of collaboration.

License
-------

This project is Copyright (c) The SunPy Developers and licensed under
the terms of the BSD 3-Clause license. This package is based upon
the `Openastronomy packaging guide <https://github.com/OpenAstronomy/packaging-guide>`_
which is licensed under the BSD 3-clause licence. See the licenses folder for
more information.

Contributing
------------

We love contributions! sunkit-instruments is open source,
built on open source, and we'd love to have you hang out in our community.

**Imposter syndrome disclaimer**: We want your help. No, really.

There may be a little voice inside your head that is telling you that you're not
ready to be an open source contributor; that your skills aren't nearly good
enough to contribute. What could you possibly offer a project like this one?

We assure you - the little voice in your head is wrong. If you can write code at
all, you can contribute code to open source. Contributing to open source
projects is a fantastic way to advance one's coding skills. Writing perfect code
isn't the measure of a good developer (that would disqualify all of us!); it's
trying to create something, making mistakes, and learning from those
mistakes. That's how we all improve, and we are happy to help others learn.

Being an open source contributor doesn't just mean writing code, either. You can
help out by writing documentation, tests, or even giving feedback about the
project (and yes - that includes giving feedback about the contribution
process). Some of these contributions may be the most valuable to the project as
a whole, because you're coming to the project with fresh eyes, so you can see
the errors and assumptions that seasoned contributors have glossed over.

Note: This disclaimer was originally written by
`Adrienne Lowe <https://github.com/adriennefriend>`_ for a
`PyCon talk <https://www.youtube.com/watch?v=6Uj746j9Heo>`_, and was adapted by
sunkit-instruments based on its use in the README file for the
`MetPy project <https://github.com/Unidata/MetPy>`_.
