%Code that shifts localizations in a slice based on the fiducial bead position 
% of the initial or previous slice (which is manually entered). This code
% also shifts targets to correct for drift between targets within a slice.

%This code is designed to take already drift corrected data for a given
%target within a slice and shift it accordingly based on the fiducial bead positions of
%previous slices and previous targets within the slice. The data inputed
%into this code should correspond to the localizations that correspond to
%only one target within one slice. 

%For this code to work, the csv should have frames as column 1, xnm as column 2,
%ynm as column 3, znm as column 4, locprec as column 5, photons as column
%6, background as column 7, locprecZ as column 8, and probe as column 9

%NS, soTILT3D 

%Import Data
clc;
clear;
[file,pathLocalizations] = uigetfile("MultiSelect","on", {'.csv'}, 'Select csv cont');
localizations_slice = readtable([pathLocalizations,file]);
localizations_array = table2array(localizations_slice);
%% Tabulate overall shifts between slices

%Enter initial positions of fiducial bead you want to shift to for LB1
%(which was the first target taken)
bead_x0_LB1 = 30085.25; %initial FB position of LB1 of 0 slice for x  
bead_y0_LB1 = 44579.34; %initial FB position of LB1 of 0 slice for y
bead_z0_LB1 = 1245.7; %initial FB position of LB1 of 0 slice for z

%Enter initial positions of fiducial bead for current slice whose localizations you want to
%shift for LB1 (which was the first target taken)
bead_x1_LB1 = 30215.74; %FB position of LB1 of current slice for x  
bead_y1_LB1 = 44838.95; %FB position of LB1 of current slice for y
bead_z1_LB1 = -5.1465; %FB position of LB1 of current slice for z  

%Tabulate shift between current slice and slice you want to shift to
x_shift = bead_x0_LB1 - bead_x1_LB1; 
y_shift = bead_y0_LB1 - bead_y1_LB1; 
z_shift = bead_z0_LB1 - bead_z1_LB1; 

%% Tabulate shifts between targets within the slice

%Enter initial positions of fiducial bead for current slice for all three
%targets 
bead_x_0_LB1 = 17649.97; %initial FB position of LB1 within the slice in x
bead_y_0_LB1 = 8080.428; %initial FB position of LB1 within the slice in y
bead_z_0_LB1 = -733.43; %initial FB position of LB1 within the slice in z
bead_x_0_LAP2 = 17665.5; %initial FB position of LAP2 within the slice in x
bead_y_0_LAP2 = 8320.53; %initial FB position of LAP2 within the slice in y
bead_z_0_LAP2 = -732.57; %initial FB position of LAP2 within the slice in z
bead_x_0_LAC = 17492.63; %initial FB position of LAC within the slice in x
bead_y_0_LAC = 8573.462; %initial FB position of LAC within the slice in y
bead_z_0_LAC = -732.09; %initial FB position of LAC within the slice in z

%Tabulate shift between second (LAP2) and third (LAC) target from first
%target (LB1) within the slice
%LAP2
x_shift_LAP2 = bead_x_0_LB1 - bead_x_0_LAP2;
y_shift_LAP2 = bead_y_0_LB1 - bead_y_0_LAP2;
z_shift_LAP2 = bead_z_0_LB1 - bead_z_0_LAP2;
%LAC
x_shift_LAC = bead_x_0_LB1 - bead_x_0_LAP2;
y_shift_LAC = bead_y_0_LB1 - bead_y_0_LAP2;
z_shift_LAC = bead_z_0_LB1 - bead_z_0_LAP2;

%% Tabulate total shift (note that x axis is inverted between Easy DHPSF and DECODE) 
z_shift_slice = 1000; %this value accounts for the z shift between slices (1000 for +1 slice, 2000 for +2 slice, etc.) 

%LB1 total shift
LB1_totalshift_x = -x_shift;
LB1_totalshift_y = y_shift;
LB1_totalshift_z = z_shift + z_shift_slice; 

%LAP2 total shift
LAP2_totalshift_x = -x_shift + (-x_shift_LAP2);
LAP2_totalshift_y = y_shift + y_shift_LAP2;
LAP2_totalshift_z = z_shift + z_shift_LAP2 + z_shift_slice; 

%LAC total shift
LAC_totalshift_x = -x_shift + (-x_shift_LAC);
LAC_totalshift_y = y_shift + y_shift_LAC;
LAC_totalshift_z = z_shift + z_shift_LAC + z_shift_slice; 

%% Apply shift to corresponding imported csv (comment or uncomment code as needed depending on target) 

%Remove all z values greater than 1 um to remove aberrations from
%localizations beyond than PSF z range
localizations_array = localizations_array(localizations_array(:,4) <=1000, :);

%LB1 (comment this block if shifting LAP2 or LAC) 
localizations_array(:,2) = localizations_array(:,2) + LB1_totalshift_x; 
localizations_array(:,3) = localizations_array(:,3) + LB1_totalshift_y; 
localizations_array(:,4) = localizations_array(:,4) + LB1_totalshift_z; 

% %LAP2 (comment this block if shifting LB1 or LAC) 
% localizations_array(:,2) = localizations_array(:,2) + LAP2_totalshift_x; 
% localizations_array(:,3) = localizations_array(:,3) + LAP2_totalshift_y; 
% localizations_array(:,4) = localizations_array(:,4) + LAP2_totalshift_z; 

% %LAC (comment this block if shifting LAP2 or LAC) 
% localizations_array(:,2) = localizations_array(:,2) + LAC_totalshift_x; 
% localizations_array(:,3) = localizations_array(:,3) + LAC_totalshift_y; 
% localizations_array(:,4) = localizations_array(:,4) + LAC_totalshift_z; 


%% Export new shifted csv

%Make table for new csv (export this table if you want to save the new
%.csv)
shifted_localizations = array2table(localizations_array);
shifted_localizations.Properties.VariableNames{1} = 'frames';
shifted_localizations.Properties.VariableNames{2} = 'xnm';
shifted_localizations.Properties.VariableNames{3} = 'ynm';
shifted_localizations.Properties.VariableNames{4} = 'znm';
shifted_localizations.Properties.VariableNames{5} = 'locprec';
shifted_localizations.Properties.VariableNames{6} = 'photons';
shifted_localizations.Properties.VariableNames{7} = 'background';
shifted_localizations.Properties.VariableNames{8} = 'locprecZ';
shifted_localizations.Properties.VariableNames{9} = 'probe';
writetable(shifted_localizations,'C:\Users\nms6\Desktop\20230713_3D LB1,LAC,LAP2 ANALYSIS\Analysis Codes and Procedures_NATCOMM\SliceStitching\Stitching_FBposition\DC,shifted,LB1,+1umslice.csv'); %be sure to change path to desired directory and name of file
