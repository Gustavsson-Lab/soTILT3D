soTILT3D
For slice stitching: 

1. Localizations in each slice were manually moved based on fiducial bead positions
See code in folder "Stitching_FBposition" and corresponding "README.txt" file

2. Localizations in each slice were better aligned in z using a cross correlation code
See code in folder "Cross_correlation_z" and corresponding "README.txt" file

3. The tops of each slice were faded based on a localization precision graded to remove harsh stitching lines. 
For top slices, this fading was performed on the bottom instead of the top. 
See codes in folder "Fading_locs_avoid_stitching_lines" and corresponding "README.txt" file

All codes here should finish running in a few minutes. Cross correlation code may take 10+ minutes.