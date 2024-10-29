soTILT3D
For line scan analysis on a 500 nm slice of the nucleus for all three nuclear targets

For analysis: 
1. Line scans were taken across all three targets and fit to a Gaussian in ImageJ. Gaussian fits were saved as csvs. 
- to do this, the same 500 nm slice of all three targets was imported in ImageJ using ThunderSTORM. Line scans were
performed on each of the three targets. 
2. Linescan_graphs.m was used to save a seperate csv for x and y for each line scan for each target. 
Three sets of all line scan datasets were saved: one that was zeroed to LAC, one that was zeroed to LAP2, 
and one that was zeroed to LB1 
3. LaminACtoLAP2.m, LaminB1toLaminAC.m, LaminB1toLAP2.m all: 
- import all the x and y csv's for the Gaussian fits for all 25 line scans for their corresponding 2 targets
- takes the average distribution for each target across all line scans
- fits the averages to Gaussian distributions and plots them
- determines the peaks for all the datasets for the two targets and determines the distances between them with 
corresponding statistics and box and whisker plots

For the Linescan_graphs.m code, the .csv exported from ImageJ for each of the three targets corresponding to a certain line scan
needs to first be Imported in Matlab using the "Import Data" button in the Home tab. Then, the code can be run. Lines 30, 31, and
32 can be edited accordingly to zero the data to the desired target. 

For the LaminACtoLAP2.m, LaminB1toLaminAC.m, and LaminB1toLAP2.m codes, the paths to pull line scan data from in lines
25, 46, 67, and 88 can all be changed to the directory that has the line scans.The prefixes can also be changed in 
lines 28, 49, 70, and 91 if desired. 

The data in the folders "ALL LINE SCANS - ZEROED TO LAC", "ALL LINE SCANS - ZEROED TO LAP2",  and "ALL LINE SCANS - ZEROED TO LB1"
was used in conjunction with the codes listed to generate statistics and plots corresponding to line scans. 