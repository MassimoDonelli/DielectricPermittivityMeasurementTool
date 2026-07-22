# About the project
EpsMeasurement.py is a tool that permits measuring the dielectric permittivity of liquid and solid materials by using an open-ended coaxial probe and a network vector analyzer; in particular, the low-cost NanoVNA V1 has been used in this example. As a probe, we used an open-ended section of semi-rigid coaxial cable equipped with a subminiature type A connector SMA.
In the source directory it is a calibration file Calibration_M.txt. In the EXAMPLE dir there are some examples of measurements with this tool, while in the CoaxProbeDetails you can find a couple of photos of the probe. 
It works for a minimum frequency of 500MHz the eps estimation has been obtained by following the guidelines reported in the reference paper:

"Dielectric Permittivity Measurement Using Open-Ended Coaxial Probe - Modeling and Simulation Based on the Simple Capacitive Model properties of minerals at microwave heating frequencies using an open-ended coaxial line," by Antonio Saloric, and Andela Matkovic, Sensor DOI 10.3390/s22166024



The code is written in python3.


The tool permits calibrating an open probe or making a measurement. To measure, you need a calibration file previously created by means of the same tool.

For the calibration, you need to perform an open-end measurement. A short-end measurement; we did it by immersing the head of the probe in mercury to obtain a good short. Then, finally, a water sample at a temperature of 25 °C.

When performing the measure especially with liquid, pay attention to air bubbles that can strongly influence the measurement.

After the measurement, you can save a graphical representation or a txt file.


Usage: python3 Source/EpsMeasurement.py 


License:
See LICENSE.txt for more information.

Contact: Massimo Donelli - massimo.donelli@unitn.it 

Project link: 