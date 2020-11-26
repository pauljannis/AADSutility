# AADS utility
Utility for absorbance-activated droplet sorting (AADS)


## Live sorting with Arduino

Code to run an Arudino Due microcontroller for AADS as updated in [Zurek et al.](https://doi.org/10.1039/D0LC00830C "Growth amplification in ultrahigh-throughput microdroplet screening increases sensitivity of clonal enzyme assays and minimizes phenotypic variation") (2020) _Lab on a Chip_.

To run this script, connect the detector to the input at pin A0 and the sorting electronic to pin D13 of the arduino.

## Offline analysis with Python

Code to analyse recorded absorabnce traces offline, e.g. to filter out fusions and create pretty histograms, can be run with the supplied python scripts.


