%Code that filters out localization precision based on a gradient depending
%on slice size and localization precision. Filtering will happen based on a
%loop that will filter for -1 nm in localization precision for every slice
%of a specified height until it reaches the max localization precision that
%the user sets.
%ADAPTION OF fading_by_LP_withloop.m BUT FOR THE BOTTOM OF THE TOP SLICE
%INSTEAD OF THE TOP OF THE BOTTOM SLICE
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
FB = 0; % set this value to a z value that is higher than your FB but lower than your data so the FB can be filterd out
A = A(A(:,4) > FB, :); % removes any fiducial beads by removing anything with a z value lower than the value set for FB
A = A(A(:,4) < max(B(:,4)), :); %ensures bottom slice isn't taller than top slice
B = B(B(:,4) > min(A(:,4)), :); %ensures bottom of top slice isn't further down than bottom of bottom slice
B_nm1 = [];
B_nm = [];

%Looping to remove based on LP on a gradient
SliceSize = 30; % determines height of each slice that is taken as loop is iterated
LP = 30; % determines LP that loop filters for as it iterates 
for i = 1:LP
    B_nm1 = B(B(:,4) <= (min(B(:,4) + (i*SliceSize))), :);
    B_nm1 = B_nm1(B_nm1(:,4) >= (min(B(:,4) + ((i-1)*SliceSize))), :);
    B_nm1 = B_nm1(B_nm1(:,5) <= i+1, :);
    B_nm = [B_nm; B_nm1];
end


%Create new array with top missing of original and add new array made in
%loop with new filtered top
B_data = B(B(:,4) > (min(B(:,4) + (LP*SliceSize))), :);
concat_newslice = [B_data; B_nm];
concat_topbottom = [A; concat_newslice];
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

%Make table for top slice (export this table if you want to save the new
%.csv)
newtopslice = array2table(concat_newslice);
newtopslice.Properties.VariableNames{1} = 'frames';
newtopslice.Properties.VariableNames{2} = 'xnm';
newtopslice.Properties.VariableNames{3} = 'ynm';
newtopslice.Properties.VariableNames{4} = 'znm';
newtopslice.Properties.VariableNames{5} = 'locprec';
newtopslice.Properties.VariableNames{6} = 'photons';
newtopslice.Properties.VariableNames{7} = 'background';
newtopslice.Properties.VariableNames{8} = 'locprecZ';
newtopslice.Properties.VariableNames{9} = 'probe';
