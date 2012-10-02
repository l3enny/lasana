# Laser Absorption Spectroscopy ANalysis

## Introduction
This is a set of utilities meant to simplify the problem of analyzing laser absorption spectra. Its present development is determined by its use in measuring the absolute, line-integrated, densities of helium metastable atoms. It should also be able to analyze the temperature (via Doppler broadening), gas velocity (via Doppler shifts), and properly account for pressure broadening. While this is only a relatively simple profile analysis, it should be sufficiently accurate for basic diagnostics.

## Dependencies
I'll try to keep these to a minimum. It's unlikely that this script uses any cutting edge features, so the version requirements are probably flexible. So far, there's:

* Numpy 1.6+
* Scipy 0.11+
* Matplotlib 1.1+

## Getting Started
Heavy lifting is handled by the `analyze` module and its submodules. `parse` and `preprocess` are both user-specific files which should be edited to suit your needs. The first is the interpreter which grabs the data files and converts them into Numpy arrays. The second applies any necessary preprocessing corrections (background subtraction, baseline shifts, that sort of thing). 