CODEX data overview
====================

CODEX is an imaging mission and most data from the mission are images.  Because the corona and solar wind are very faint
compared to the various backgrounds in the data, high photometric precision and many steps of processing are required.  This
means finding CODEX data can be complicated.

CODEX is designed to collect images of the solar corona in
polarized light within the wavelength range spanning from 385 to
440 nm to measure the coronal electron density, temperature, and
speed between 3 and 8 solar radii. CODEX is a polarimeter, so image data arrive in polarized form, consisting
of four-directional (0, 45, 90, and 135°) pixels. The CMOS has 4096 × 3000 pixels of 3.45μm x 3.45μm. CODEX
utilizes this detector to simultaneously image each of the four polarizer angles in
a single observation, therefore eliminating any evolution of the scene between polarization observations.


The 4 levels of processing are:

Level 0
-------
These are data direct from the CODEX camera, assembled into FITS files with 16 bit pixels. Headers only contain information that came from instrument. The data are typically lossy compressed on board utilizing JPEG-LS; L0 images
have been decompressed into their original form, values given in raw counts (DN). This data is currently not provided publicly.

Level 1
-------
These are the FITS data, merged with additonal header information and split into the 4 polarization states. The data uses FITS extensions. Values are in raw counts (DN). Header contains all information in telemetry plus any ancillary information necessary to interpret the data, for example full World Coordinte System (WCS) compliance.

Level 2
-------
FITS files from a single detector with calibrations applied. Values are in physical units (solar brightness).


Level 3
-------
Data products are the result of combining two or more images (mosaics, movies, Carrington maps, etc.) or Derived quantities (pB, electron densities, temperature, and velocity, etc.). May or may not be physical units.

