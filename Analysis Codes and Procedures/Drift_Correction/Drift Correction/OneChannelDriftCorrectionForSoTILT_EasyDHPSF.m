%% Drift correction code for soTILT - Easy DHPSF bead and DECODE data
% Code for drift correction of 1 channel, 3 dimenesional data from DECODE and a
% transformed bead from Easy DHPSF. Edited 3/9/23 by JL for soTILT project to be Easy DHPSF compatibale , 
% adapting material from TransformSMACMData_v6.m
% Input variables from SMAP in csv are xnm, ynm, znm, frame, locprecnm,
% photon, bg, locprecznm ( in this order )
% *** indicates that the code will need to be edited here between samples
clear all
clc
close all

%% Load data from SMAP
[file,path] = uigetfile("MultiSelect","on", {'.csv'}, 'Select all SMACM data from SMAP');
% Reads the tables and updates their frames for all files 
SMACM_data = [];
cell_test = iscell(file);
if cell_test == 1
for i = 1:length(file)
    newTable = readtable([path, file{i}]);
    if i ==1
    else
        lastFrame = SMACM_data(end,4);
        lastFrame = table2array(lastFrame);
        newTableFrames = table2array(newTable(:,4));
        newTableFrames = newTableFrames+ lastFrame;
        newTableFrames = array2table(newTableFrames);
        newTable(:,4) = newTableFrames;
    end
    SMACM_data = [SMACM_data; newTable];
    clear newTable
end
else
SMACM_data = readtable([path, file]);
end


%% Fiducial bead in another csv
% If the FB data is in another csv, uncomment this section
[fileFB, pathFB] = uigetfile("MultiSelect","on", {'.csv'}, 'Select fiducial bead data from SMAP');
fiducial_localizations = readtable([pathFB,fileFB]);
SMACM_data = table2array(SMACM_data);
fiducial_localizations = table2array(fiducial_localizations);
Drift_correction_data = [SMACM_data(:,4) SMACM_data(:,1) SMACM_data(:,2) SMACM_data(:,3), SMACM_data(:,5) SMACM_data(:,6) SMACM_data(:,7) SMACM_data(:,8); fiducial_localizations(:,1) fiducial_localizations(:,3) fiducial_localizations(:,4) fiducial_localizations(:,5) fiducial_localizations(:,2) fiducial_localizations(:,6) fiducial_localizations(:,7) fiducial_localizations(:,2)];
SMACM_data = array2table(SMACM_data);
%SMACM_data = [SMACM_data;fiducial_localizations];
% SMACM_data = table2array(SMACM_data);
% SMACM_data = sortrows(SMACM_data);
% SMACM_data = array2table(SMACM_data);
%SMACM_data_variablenames = {'frame', 'x_pix','y_pix','z_nm','photons','background','crlb_x','crlb_y','crlb_z','crlb_photons','crlb_background','logLikelyhood','x_nm','y_nm', 'crlb_xnm','crlb_ynm'};
% SMACM_data.Properties.VariableNames = SMACM_data_variablenames;

%% Track and isolate the fiducial bead 
% Approximate initial bead localizations
xBead = 7;%px
yBead = 7;%px

xBead_nm = xBead*159;
yBead_nm = yBead*159;

xnm = Drift_correction_data(:,2);
ynm = Drift_correction_data(:,3);
znm = Drift_correction_data(:,4);
frames = Drift_correction_data(:,1);

% Plot data to find bead 
figure (1)
scatter(xnm,ynm,[],frames,'filled');
set(gca,'YDir','reverse')
% Coloring
h2 = colorbar;
colormap jet
ylabel(h2,'frames [#]')
% bounds of the plot
xlim([(xBead_nm-2000) (xBead_nm+2000)])
ylim([(yBead_nm-2000) (yBead_nm+2000)])
title('Circle Bead Localizatons')

% Circle bead 
hpolygon = drawfreehand();
[in,on] = inpolygon(xnm,ynm,hpolygon.Position(:,1),hpolygon.Position(:,2));

% Isolate bead localizations while keeping index 
xBead = zeros(length(xnm),2);
yBead = zeros(length(ynm),2);
zBead = zeros(length(znm),2);

for i = 1:length(xnm)
    xBead(i,:) = [xnm(i)*in(i), frames(i)];
    yBead(i,:) = [ynm(i)*in(i), frames(i)];
    zBead(i,:) = [znm(i)*in(i), frames(i)];
end
% Plot unfiltered bead 
figure (2)
scatter(xBead(:,1), yBead(:,1));
xlim([(xBead_nm-2000) (xBead_nm+2000)])
ylim([(yBead_nm-2000) (yBead_nm+2000)])

xBead = xBead(all(xBead,2),:);
yBead = yBead(all(yBead,2),:);
zBead = zBead(all(zBead,2),:);
%% Filter FB localizations 

[xBeadFiltered(:,2),xBeadFiltered(:,1)] = filterlocs(xBead(:,2),xBead(:,1),2000,1,'x','D:\EasyDHPSF & DECODE Validation\DECODE\Bead 1\Bead1 drift correction\xBeadFiltered.mat');
[yBeadFiltered(:,2),yBeadFiltered(:,1)] = filterlocs(yBead(:,2),yBead(:,1),2000,1,'y','D:\EasyDHPSF & DECODE Validation\DECODE\Bead 1\Bead1 drift correction\yBeadFiltered.mat');
[zBeadFiltered(:,2),zBeadFiltered(:,1)] = filterlocs(zBead(:,2),zBead(:,1),2000,1,'z','D:\EasyDHPSF & DECODE Validation\DECODE\Bead 1\Bead1 drift correction\zBeadFiltered.mat');

%% Track Fiducial bead 

beadInitial = [xBeadFiltered(1,1) yBeadFiltered(1,1) zBeadFiltered(1,1)];
save("Initial bead position Im1.mat", "beadInitial", '-mat');

% Cubic smoothing function 
fidTrackX = csaps(xBeadFiltered(:,2),xBeadFiltered(:,1), 0.0000000001);
figure(4)
fnplt(fidTrackX);
hold on
scatter(xBeadFiltered(:,2) , xBeadFiltered(:,1));

beadArrayX = fnval(fidTrackX, frames);

fidTrackY = csaps(yBeadFiltered(:,2),yBeadFiltered(:,1), 0.0000000001);
figure(5)
fnplt(fidTrackY);
hold on
scatter(yBeadFiltered(:,2) , yBeadFiltered(:,1));

beadArrayY = fnval(fidTrackY, frames);

fidTrackZ = csaps(zBeadFiltered(:,2),zBeadFiltered(:,1), 0.0000000001);
figure(6)
fnplt(fidTrackZ);
hold on
scatter(zBeadFiltered(:,2) , zBeadFiltered(:,1));

beadArrayZ = fnval(fidTrackZ, frames);

% Deviation in the bead from the inital bead position 
 for i = 1:length(frames)
     devBeadX(i) = beadArrayX(i) - beadArrayX(1);
     devBeadY(i) = beadArrayY(i) - beadArrayY(1);
     devBeadZ(i) = beadArrayZ(i) - beadArrayZ(1);
 end
 devBeadX = devBeadX';
 devBeadY = devBeadY';
 devBeadZ = devBeadZ';

 %% Apply fiducual shift to all data 

 for i = 1:length(frames)
     xnmDriftCorrected(i) = xnm(i)+devBeadX(i);
     ynmDriftCorrected(i) = ynm(i)-devBeadY(i);
     znmDriftCorrected(i) = znm(i)-devBeadZ(i);
 end

xnmDriftCorrected_tab = array2table((xnmDriftCorrected'));
xnmDriftCorrected_tab = renamevars(xnmDriftCorrected_tab, "Var1", 'xnm');
ynmDriftCorrected_tab = array2table((ynmDriftCorrected'));
ynmDriftCorrected_tab = renamevars(ynmDriftCorrected_tab, "Var1", 'ynm');
znmDriftCorrected_tab = array2table((znmDriftCorrected'));
znmDriftCorrected_tab = renamevars(znmDriftCorrected_tab, "Var1", 'znm');

probe = zeros(length(frames), 1);
probe = array2table(probe);
frames = array2table(frames);
Drift_correction_data = array2table(Drift_correction_data);
% Format final csv file 
Formatted_DriftCorrectedSMACMData = [ frames xnmDriftCorrected_tab ynmDriftCorrected_tab znmDriftCorrected_tab  Drift_correction_data(:,5) Drift_correction_data(:,6) Drift_correction_data(:,7) Drift_correction_data(:,8) probe];
Formatted_DriftCorrectedSMACMData.Properties.VariableNames = {'frames', 'xnm','ynm','znm','locprec','photons','background','locprec z','probe'};

writetable(Formatted_DriftCorrectedSMACMData, 'D:\EasyDHPSF & DECODE Validation\DECODE\Bead 1\Bead1 drift correction\Bead1_driftcorrected.csv');

