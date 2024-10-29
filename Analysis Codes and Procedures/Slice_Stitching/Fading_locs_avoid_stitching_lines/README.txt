Use fading_by_LP_withloop_fortop.m if you want to filter the top of hte bottom slice. 
You can vary depth and what localization precision to filter for. This will filter the
data based on a gradient where the LP filtering gets more agressive towards the topmost 
part of the slice. 

Use fading_by_LP_withloop_forbottom.m if you want to filter the bottom of the top slice. 
This should only be relevant for the topmost slices of the nucleus. You can vary depth and
what localization precision to filter for. This will filter the data based on a gradient 
where the LP filtering gets more aggressive towards the bottommost part of the slice. 

For this code to work, the csv should have frames as column 1, xnm as column 2, ynm as column 3,
znm as column 4, locprec as column 5, photons as column 6, background as column 7, locprecZ as 
column 8, and probe as column 9
If using the Shiftinglocs_basedon_FBpositions.m code prior, this data should already be in the correct format

In this code, you: 

1. Import the data. Select the correct csv file when prompted for top and bottom slices. 

2. Loop through the localizations to remove localizations based on a localization
precision gradient. The value specified in line 31 is the SliceSize and the value in line
32 is the localization precision that the loop filters for as it iterates. The code will loop 
through and remove localizations based on decreasing localization precision, with the lowest localization
precision being at the top and the highest localization precision being towards the bottom. The depth
through which it removes localizations for each localization precision value is determined by the SliceSize. 
This value can be modified if you want more gradient filtering with a greater or lesser depth. The value for LP
is typically chosen to be the same value that all the data is filtered for/the max localization precision in the dataset. 

3. Create two figures. One figure is the YZ profile of the two slices before and the YZ profile of the two slices after 
gradient filtering. This can be done to assess how well the code has removed the harsh stitching lines. 

4. Make a new csv. 
If satisfied with the gradient filtering, the new csv can be exported using the "writetable" function in Matlab. 

To use this code, the following lines may need to be manually modified: 
31, 32 

The fading_by_LP_withloop_fortop.m code can be tested with the sample input data provided for the top and bottom slices. 
The sample output data corresponding to the input data is also provided. 
The fading_by_LP_withloop_forbottom.m is just an analog of the fortop version, but filters the bottom of the top slice
instead of the top of the bottom slice. 