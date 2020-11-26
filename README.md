# AADS utility
Utility for absorbance-activated droplet sorting (AADS)


## Live sorting with Arduino

Code to run an Arudino Due microcontroller for AADS as updated in [Zurek et al.](https://doi.org/10.1039/D0LC00830C "Growth amplification in ultrahigh-throughput microdroplet screening increases sensitivity of clonal enzyme assays and minimizes phenotypic variation") (2020) _Lab on a Chip_.

To run the ArduinoAADS.ino script, connect the detector to the input at pin A0 and the sorting electronic to pin D13 of the arduino.
Now set the sorting variables:
- thresh: Voltage threshold for peak detection. Set this ~0.5V below the baseline.
The 2D sorting gate is determined by the peak height and peak width. Droplets in the range sortVH-sortVL and peakWH-peakWL will thus be sorted.
- sortVH: High sorting voltage threshold. Peaks lower than this value will be sorted.
- sortVL: Low sorting voltage threshold. Peaks higher than this value will be sorted.
- peakWH: High size threshold. Smaller droplets than this value will be sorted.
- peakWL: Low size threshold. Droplets larger than this value will be sorted.


## Offline analysis with Python

Code to analyse recorded absorbance traces offline, e.g. to filter out fusions and create pretty histograms, can be run with the supplied python scripts. Running these scripts requires python 3 with some libraries (numpy, matplotlib, seaborn).

### Stage 1: AADS_detect.py
Run this script to extract peak information, such as peak height and width, from a raw recording of detection signal. The detection signal needs to be in the format of voltage values with timestamps in seconds, as given in example-trace.txt. Run the peak detection with `python AADS_detect.py midpoint example-trace.txt -o` or consult the help `python AADS_detect.py -h`. 

The output will be a visualization of the trace, as in the example below, and a file with the extracted peaks.



