Introduction
=============

What is CODEX?
--------------
Understanding solar wind sources and acceleration mechanisms is an overarching solar physics goal. Current models are highly under-constrained due to the limitations of the existing data, particularly in the ~3-10 Rs range. COronal Diagnostic EXperiment (CODEX) is designed to deliver the first global, comprehensive data sets that will impose crucial constraints and answer targeted essential questions, including: Are there signatures of hot plasma released into the solar wind from previously closed fields? What are the velocities and temperatures of the density structures that are observed so ubiquitously within streamers and coronal holes? To provide these crucial measurements, NASA’s Goddard Space Flight Center, in collaboration with the Korea Astronomy and Space Science Institute, and the Italian National Institute for Astrophysics have developed a next-generation coronagraph hosted on the International Space Station, launched November 4, 2024. This imaging solar coronagraph uses multiple narrow passband filters to obtain simultaneous measurements of electron density, temperature, and velocity of the nascent solar wind. CODEX provides comprehensive data sets that test theories of solar wind formation and provide crucial constraints on predictive solar wind models, taking multiple daily measurements in this critical solar wind formation region. See https://science.nasa.gov/mission/codex/ for more details.

How do I get data?
-----------------

You can learn more at `this page <data/index.html>`_.

Where does `data-pipeline` fit in?
------------------------------
``data-pipeline`` is the data reduction pipeline code for the CODEX mission. The pipeline, as shown below,
consists of several segments of processing.

These segments are the following:

- Raw to *Level 0*: converts raw coronagraph data to FITS images
- Level 0 to *Level 1*: Updates FITS headers to include all necessary informaiton for science analysis, splits data into 4 polarization states.
- Level 1 to *Level 2*: Provides calibrated, polarized images
- Level 2 to Level 3: higher level products, combined images.

CODEX and Python
----------------

The CODEX framework is built using Python - an object-oriented language with a large user / code base in astronomy and solar physics.
The pipeline and tools for querying / loading PUNCH data use the Python language, along with the SunPy and Astropy software libraries.
A number of useful tutorials exist online...

In addition to scripts and modules, Python notebooks provide a great way to execute and document a sequence of code cells, with visualizations directly in-line.
It's a sort of analogue of the classic research notebook.
