CODEX data structure
====================

Concept
-------
CODEX data is structured using the ??? framework, which bundles an N-dimensional data array with a corresponding World Coordinate System (WCS) describing the spatial or spectral coordinates of the data itself.
This allows for better integration with AstroPy software libraries, including coordinate and data reprojection or resampling.
For primary Level 2 CODEX data, connecting data with coordinates is critical given the relatively large field of view of the combined virtual CODEX observatory.

Uncertainty
-----------

The uncertainty is stored within the ???