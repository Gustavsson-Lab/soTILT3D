This code takes two slices (a top and bottom slice) and performs a z cross correlation. This code is used after the slices are manually shifted
using a fiducial bead. This code helps correct for any residual mismatch between slices in z after FB shifting. It does so by optimizing the
overlap between the top 200 nm of the bottom slice and the bottom 200 nm of the top slice. 

To run this code, the Parallel Computing Toolbox must be installed on Matlab and the subfunction "myfuncc.m" must be in the same folder as the code

1. csv files corresponding to the bottom and top slice that you would like to cross correlate must be imported into matlab
You can do this by clicking Home --> Import Data to import the correct csv's of the bottom and top slice that you wish to cross correlate into the workspace
Turns each csv into an array in the commandwindow by using the "table2array" function in matlab

2. The distance of overlap between the two slices that the cross correlation is tabulated over can be modified in line 30 if desired (default set to 200 nm) 

2. After running code, the Z correction value will appear in the workspace as a single value with the variable name "Zcorrection". Figures will also pop up
showing the Before and After Z correction so the user can access the performance of the cross correlation by eye. A shift in z corresponding to this value should then be 
applied manually the z values in the csv. 

To use this code, the following lines may need to be manually modified: 
30

**NOTE: This code does not give an output csv. This code gives a value corresponding to a z shift that must be manually performed on the data.
In some cases, the data in the slices may need to be "cropped" (ex: taking only the bottom 500 nm of the top slice and the top 500 nm of the top slice rather than the entire 2 micron range)
to optimize overlap.  