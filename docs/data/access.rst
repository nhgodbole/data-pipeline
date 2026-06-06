Accessing CODEX Data
====================

Recommended Level 1 Products
-----------------------------
For most science use cases, the recommended starting point is the Level 1 science data products...


Downloading Data
----------------
Data output from the CODEX data processing pipeline are stored and accessible through the Solar Data Analysis Center (SDAC)
- a portal for hosting through tools such as the Virtual Solar Observatory (VSO).
From here CODEX data products can be queried and requested for download using metadata within the data products.

If that example is not working properly, you can also pull data directly from the SDAC using ``wget``.

.. code-block:: bash

    wget -r -l1 --no-parent --no-directories -A "codex_l1_20250521_*.fits" -R "*.html*,index*,*tmp*" https://umbra.nascom.nasa.gov/codex/2025/05/21/

The above example would pull data for the L1 products on 2025-05-21.
Change the path and date according to what product you wish to download.

CODEX data are also accessible using the ??? tool,
where it can be quickly pulled.

Reading Data
------------
Standard CODEX data is stored as a standards-compliant FITS file, which bundles the primary data along with secondary data and metadata fully describing the observation.
Each file is named with a convention that uniquely identifies the product - a sample being 'codex_l1_20250521_000001_5_4.fits' - where l1 defines the data level,
20250521_000001 is a timestamp in the format yyyymmdd_hhmmss, and _5, _4 are the positions of the filter wheels.

For most end-users the primary data of interest are L1.

These data are compatible with standard astropy FITS libraries, and can be read in as following the example,

.. code-block:: python

    filename = 'example_data/codex_l1_20250521_000001_5_4.fits.fits'

    with fits.open(filename) as hdul:
        data = hdul[1].data
        header = hdul[1].header
        uncertainty = hdul[2].data

Data Projections
----------------
The CODEX instrument extends its field of view out to around ???-degrees from the Sun,
creating a meshed virtual observatory extending to a diameter of nearly ??? solar radii.
The wide nature of this field of view requires attention to the data projection being used for these data.

Each data contains a set of World Coordinate System (WCS) parameters that describe the coordinates of the data,
in both a helioprojective and celestial frame.
