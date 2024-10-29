%% 2D Transformation for FB - EASY DHPSF
%% load data 
%Drift_Corrected_data;
[file,path] = uigetfile("MultiSelect","on", {'.csv'}, 'Select FB data from Easy DHPSF');
FB_data = [];

% Sorts FB localizations in order of acquisition and compiles them into one
% array with accurate frame numbers
cell_test = iscell(file);
if cell_test == 1
    for i = 1:length(file)
    newTable = readtable([path, file{i}]);
    if i ==1
    else
        lastFrame = FB_data(end,1);
        lastFrame = table2array(lastFrame);
        newTableFrames = table2array(newTable(:,1));
        newTableFrames = newTableFrames+ lastFrame;
        newTableFrames = array2table(newTableFrames);
        newTable(:,1) = newTableFrames;
    end
    FB_data = [FB_data; newTable];
    clear newTable
    end
else
end
FB_data = readtable([path, file]);
FB_data = table2array(FB_data);
% take out the NaN data points 
FB_data = [FB_data(:,1) FB_data(:,2) FB_data(:,3) FB_data(:,4) FB_data(:,5) FB_data(:,9) FB_data(:,10)];

Drift_Corrected_data_variablenames = {'frame','molecule' ,'x_nm','y_nm','z_nm','photons','background'};
% Scatterplot drift corrected data
figure (1)
scatter(FB_data(:,3), FB_data(:,4),'r.');
% get tform.mat
tform=uigetfile('.mat', 'Open tform.mat file');
load(tform);

x_px = FB_data(:,3)./159; 
y_px = FB_data(:,4)./159;

[transformedx,transformedy] = transformPointsForward(tform, x_px,y_px);
figure (2)
scatter(transformedx,transformedy,'g.');
title('Drift Corrected and transformed 2D');
hold off
transformedx = transformedx.*159;
transformedy = transformedy.*159;
FB_Transformed = FB_data;
FB_Transformed(:,3) = transformedx;
FB_Transformed(:,4) = transformedy;
%folderpath = 'C:\Users\LabUser\Desktop\20230312_DECODETEST\TESTING_EASYDHPSFBEAD_DATADECODE';
FB_Transformed = array2table(FB_Transformed);
FB_Transformed.Properties.VariableNames = Drift_Corrected_data_variablenames;
writetable(FB_Transformed, 'C:\Users\gpg2\Documents\20230630 Microtubules, Mitochondria, Lap2\0.08 nM Lap2\Slice 0\FB\FB transformed_Red to Green_Lap2_0um.csv');

