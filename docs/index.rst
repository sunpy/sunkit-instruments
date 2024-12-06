********************************
sunkit-instruments Documentation
********************************

A package for instrument-specific data structures and processing in the SunPy ecosystem.

Installation
============

For detailed installation instructions, see the `installation guide`_ in the ``sunpy`` docs.
This takes you through the options for getting a virtual environment and installing ``sunpy``.
You will need to replace "sunpy" with "sunkit-instruments".

Getting Help
============

Stop by our `chat room <https://app.element.io/#/room/#sunpy:openastronomy.org>`__ if you have any questions.

Contributing
============

Help is always welcome so let us know what you like to work on, or check out the `issues page`_ for the list of known outstanding items.
If you would like to get involved, please read our `contributing guide`_, this talks about ``sunpy`` but the same is for ``sunkit-instruments``.

If you want help develop ``sunkit-instruments`` you will need to install it from GitHub.
The best way to do this is to create a new python virtual environment.
Once you have that virtual environment, you will want to fork the repo and then run::

    $ git clone https://github.com/<your_username>/sunkit-instruments.git
    $ cd sunkit-instruments
    $ pip install -e ".[dev]"

.. _installation guide: https://docs.sunpy.org/en/stable/tutorial/installation.html
.. _issues page: https://github.com/sunpy/sunkit-instruments/issues
.. _contributing guide: https://docs.sunpy.org/en/latest/dev_guide/contents/newcomers.html


``sunkit-instruments`` is organised into sub-modules for each instrument:

.. toctree::
   :maxdepth: 2

   code_ref/index
   generated/gallery/index
   topic_guide/index
   whatsnew/index


Note that the code in this package is **not** maintained by or necessarily contributed to by instrument teams.
Some instruments have individual packages for analysis in Python, including:

- `aiapy <https://aiapy.readthedocs.io/>`__
- `eispac <https://eispac.readthedocs.io/>`__
- `xrtpy <https://xrtpy.readthedocs.io/>`__
