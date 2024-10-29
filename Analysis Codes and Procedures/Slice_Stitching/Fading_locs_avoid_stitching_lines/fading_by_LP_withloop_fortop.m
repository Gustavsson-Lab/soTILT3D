%Code that filters out localization precision based on a gradient depending
%on slice size and localization precision. Filtering will happen based on a
%loop that will filter for -1 nm in localization precision for every slice
%of a specified height until it reaches the max localization precision that
%the user sets.

%**NOTE: REMOVE FB FROM DATA IF POSSIBLE BEFORE BEGINNING! YOU CAN TRY
%RUNNING THROUGH WITHOUT DOING THIS, BUT IF DATA LOOKS WONKY, CHECK TO SEE
%IF MAX VALUES ARE SKEWED BY A FIDUCIAL BEAD THAT IS VERY HIGH RELATIVE TO
%THE DATA

%August 16, 2023, soTILT3D 

%Import Data
clc;
clear;
[filebot,pathbottom] = uigetfile("MultiSelect","on", {'.csv'}, 'Select csv for bottom slice');
[filetop, pathtop] = uigetfile("MultiSelect","on", {'.csv'}, 'Select csv for top slice');
top_slice = readtable([pathtop,filetop]);
bottom_slice = readtable([pathbottom,filebot]);
top_slice = table2array(top_slice);
bot_slice = table2array(bottom_slice);

A = bot_slice; 
B = top_slice; 
A = A(A(:,4) < max(B(:,4)), :); %ensures bottom slice isn't taller than top slice
A_nm1 = [];
A_nm = [];

%Looping to remove based on LP on a gradient
SliceSize = 30; % determines height of each slice that is taken as loop is iterated
LP = 30; % determines LP that loop filters for as it iterates 
for i = 1:LP
    A_nm1 = A(A(:,4) >= (max(A(:,4) - (i*SliceSize))), :);
    A_nm1 = A_nm1(A_nm1(:,4) <= (max(A(:,4) - ((i-1)*SliceSize))), :);
    A_nm1 = A_nm1(A_nm1(:,5) <= i+1, :);
    A_nm = [A_nm; A_nm1];
end


%Create new array with top missing of original and add new array made in
%loop with new filtered top
A_data = A(A(:,4) < (max(A(:,4) - (LP*SliceSize))), :);
concat_newslice = [A_data; A_nm];
concat_topbottom = [concat_newslice; B];
concat_original = [A; B];

%Visualize and render
LP_s = num2str(LP);

figure
scatter(concat_original(:,3), concat_original(:,4), 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', 0.002);
title('YZ view for original')

figure
scatter(concat_topbottom(:,3), concat_topbottom(:,4), 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', 0.002);
title_YZ = strcat('YZ view for LP less than', LP_s, 'for gradient for top of bottom slice');
title(title_YZ)

%Make table for bottom slice (export this table if you want to save the new
%.csv)
newbottomslice = array2table(concat_newslice);
newbottomslice.Properties.VariableNames{1} = 'frames';
newbottomslice.Properties.VariableNames{2} = 'xnm';
newbottomslice.Properties.VariableNames{3} = 'ynm';
newbottomslice.Properties.VariableNames{4} = 'znm';
newbottomslice.Properties.VariableNames{5} = 'locprec';
newbottomslice.Properties.VariableNames{6} = 'photons';
newbottomslice.Properties.VariableNames{7} = 'background';
newbottomslice.Properties.VariableNames{8} = 'locprecZ';
newbottomslice.Properties.VariableNames{9} = 'probe';
