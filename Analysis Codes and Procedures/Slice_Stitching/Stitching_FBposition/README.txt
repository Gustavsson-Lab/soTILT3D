This code shifts localizations in a slice based on the fiducial bead position of 
the previous (or 0) slice, which is manually entered. This code
also shifts targets to correct for drift between targets within a slice.

This code is designed to take already drift corrected data for a given
target within a slice and shift it accordingly based on the fiducial bead positions of
previous slices and previous targets within the slice. The data inputed
into this code should correspond to the localizations that correspond to
only one target within one slice. 

For this code to work, the csv should have frames as column 1, xnm as column 2, ynm as column 3,
znm as column 4, locprec as column 5, photons as column 6, background as column 7, locprecZ as 
column 8, and probe as column 9
If using the drift correction code for EasyDHPSF bead correcting DECODE data, this data should
already be in the correct format

In this code, you: 

1. Import the data. Select the correct csv file when prompted. The code is designed to accept
the csv file corresponding to a specific target within a specific slice. 

2. Tabulate overall shift between slices 
This uses bead positions to tabulate the overall shift between slices. 
The initial FB positions of the first slice is manually put in in lines 27, 28, and 29. 
The FB positions of the current slice is manually put in in lines 33, 34, and 35. 

3. Tabulate shifts between targets within the slice
This uses bead positions to tabulate the overall shift between targets within a slice. 
The initial FB positions of the first target is manually put in in lines 46, 47, and 48. 
The intial FB positions of the second target is manually put in in lines 49, 50, and 51.
The initial FB positions of the third target is manually put in in lines 52, 53, and 54.

4. Tabulate total shift 
This tabulates the total shift based on the FB positions entered and the corresponding z shift between slices (entered in line 68). 

5. Apply shift to corresponding target
Comment or uncomment relevant blocks of code depending on which target you are shifting. 
In line 89, all localizations in z greater than 1 um are removed to remove aberrations from localizations beyond the PSF z range

6. Export new shifted csv
Turns array back into a table and renames columns with proper column names. Path and name of file exported is changed in line 120. 

To use this code, the following lines need to be manually modified: 
27, 28, 29, 33, 34, 35, 46, 47, 48, 49, 50, 51, 52, 53, 54, 68, 121

The code can be tested with the sample input data provided. The sample output data corresponding to the input data is also provided. 
