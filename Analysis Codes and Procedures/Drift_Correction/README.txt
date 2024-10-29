Drift correction code for single-molecule data taken on soTILT3D
Used for data acquired by Gabriella Gagliano and Nahima Saliba in the lab of Dr. Anna-Karin Gustavsson in 2022, and 2023.
(this code recognizes fiducial bead data in the red channel localized by EasyDHPSF (we used the 12 um DHPSF)
and single-molecule data in the green channel localized by DECODE (we used the 2 um range DHPSF), but can be 
easily adapted if you use different localization software which outputs localizations in a different format.
(e.g. using fiducial bead data from DECODE which such as was done to validate DECODE compared to EasyDHPSF.)
This code should run in under a few minutes.

What you will find in the folder:
1.) 12 um Easy-DHPSF software for fiducial bead localization
( You will need a calibration scan of a long range bead, a squence.mat file, a dark counts image, and your red channel fiducial bead data to run EasyDHPSF.
Examples of these are provided in the folder.)
2.) Directions to find DECODE installation instructions for single-molecule localization (Included is also a green channel high-density 3D single-molecule data set to be localized with DECODE.)
3.) Registration code for registering the two channels
( You will need registration images of TetraSpeck fiducial beads covering the whole FOV of both channels in order to obtain your transformation matrix to map the two channels.
These can be in 2D with the standard PSF. An example is provided in the folder as well as a completed Tform that will allow you to transform your data.)
4.) Transformation code for transforming red data onto green channel (You will need the Tform created in the previous step and your red channel fiducial bead localization data.
Examples of both are provided in the folder.)
5.) Drift correction code for drift correcting single-molecule data with fiducial bead data (You will need your DECODE single-molecule localization data and filterlocs.m to run this code.
Examples of both are provided in the folder.)

Directions:
1.) Crop your single-molecule data in the green channel and your fiducial bead data in the red channel the same way
2.) Localize your red channel fiducial bead data in EasyDHPSF
3.) Localize your green channel single-molecule data in DECODE
4.) Crop your registration images in both channels in the same way
5.) Register your two channels with the registration code. Be sure to change the path to save your tform.
6.) Perform the FB transformation using the FB localizations from EasyDHPSF and your tform.m file.
Make sure to change the path to save your transformed bead data.
7.) Perform the drift correction with drift correction code using your transformed bead data and your green channel single-molecule DECODE data.
Make sure filterlocs.m is in the same folder for this code to work. Input the xBead and yBead positions in pixels in line 54 and 55 based on their
positions in the transformed csv file. Remember to change the path to of the beadfit results in lines 103, 104, and 105 and the path of the final
drift corrected csv file in line 170 to where you want them saved.
8.) Filter final output data if desired to make the csv file less data heavy. This can also be done in Vutara during rendering. An "ID" column will
also need to be added if you wish to render your final data in Vutara.