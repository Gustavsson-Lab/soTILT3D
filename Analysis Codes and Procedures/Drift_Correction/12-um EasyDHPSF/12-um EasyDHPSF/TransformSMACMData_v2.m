% 2018-10-23 - CAB updated to allow separate fiducial data to be loaded in
% separately
% fiducial info should be stored in matrix called avgFid, which is an
% nFrames x 3 matrix

clear all 
close all 

nmPerPixelG = 151; nmPerPixelR = 160;
scrsz = get(0,'ScreenSize');
nSample = 1.33;         % index of refraction of sample
nOil = 1.518;           % index of immersion oil

%% load in the transform 
[tformFile, tformPath] = uigetfile({'*.mat';'*.*'},'Open 3D_Transform.mat');
if isequal(tformFile,0)
    error('User cancelled the program');
end

% tformFile = '3D_Transform_lwquadratic.mat';
% tformPath = 'C:\UserFiles\UsersNOTBackedUp\Temp\2014-04-09_L1-PAM-B-eYFP\FittingCode-Sept3rd-2014\';
% % 

load([tformPath tformFile]);

%% Prepare for data interpolation
x = matched_cp_reflected(:,5);
y = matched_cp_reflected(:,6);
z = matched_cp_reflected(:,7);
F_FRE = TriScatteredInterp(x,y,z,FRE_full(:,1), 'natural');
F_FRE_X = TriScatteredInterp(x,y,z,FRE_full(:,2), 'natural');
F_FRE_Y = TriScatteredInterp(x,y,z,FRE_full(:,3), 'natural');
F_FRE_Z = TriScatteredInterp(x,y,z,FRE_full(:,4), 'natural');
F_TRE = TriScatteredInterp(x,y,z,TRE_full(:,1), 'natural');
F_TRE_X = TriScatteredInterp(x,y,z,TRE_full(:,2), 'natural');
F_TRE_Y = TriScatteredInterp(x,y,z,TRE_full(:,3), 'natural');
F_TRE_Z = TriScatteredInterp(x,y,z,TRE_full(:,4), 'natural');


% Columns: 
frameCol = 1;
matchConCol = 6; 
amp1Col = 7; 
amp2Col = 8; 
xLoc1Col = 9; 
yLoc1Col = 10; 
xLoc2Col = 11; 
yLoc2Col = 12; 
sigma1Col = 13; 
sigma2Col = 14; 
bkgndMeanCol = 15; 
fitErrorCol = 16; 
goodFitFlagCol = 17; 
xCenterCol = 18; 
yCenterCol = 19; 
AngleCol = 20; 
PhotonCol = 21;
interlobeCol = 22; 
ampRatioCol = 23; 
sigmaRatioCol = 24; 
xCol = 25;
yCol = 26;
zCol = 27;
xFidCol = 28;
yFidCol = 29; 
zFidCol = 30; 

fileNum = 1;
LocFiles = {};
dataSets = [];

%% load in starting fiducial data
[fidStartFile, fidStartPath] = uigetfile({'*.mat';'*.*'},'Load in the starting fid data');

load([fidStartPath fidStartFile]);

startFidPosG = median(green.results.traj).*1e9; % convert to nm
startFidPosR = median(red.results.traj).*1e9; % convert to nm

% change coords from relative to absolute
startFidPosG_xShift = (min(green.data.cropCoordsC) + length(green.data.cropCoordsC)./2).*nmPerPixelG;
startFidPosG_yShift = (min(green.data.cropCoordsR) + length(green.data.cropCoordsR)./2).*nmPerPixelG;
startFidPosR_xShift = (min(red.data.cropCoordsC) + length(red.data.cropCoordsC)./2).*nmPerPixelR;
startFidPosR_yShift = (min(red.data.cropCoordsR) + length(red.data.cropCoordsR)./2).*nmPerPixelR;

startFidPosG(:,1:2) = startFidPosG(:,1:2) + [startFidPosG_xShift startFidPosG_yShift];
origStartFidPosG =  startFidPosG(:,3);
startFidPosG(:,3) = startFidPosG(:,3) - origStartFidPosG;
startFidPosR(:,1:2) = startFidPosR(:,1:2) + [startFidPosR_xShift startFidPosR_yShift];
origStartFidPosR =  startFidPosR(:,3);
startFidPosR(:,3) = startFidPosR(:,3) - origStartFidPosR;

% % transform R startingPos to G channel
% transformedStartPosR = transformData([startFidPosR(:,1),startFidPosR(:,2),startFidPosR(:,3)],tform);

%% load in the GREEN .mat file
% Add it to the total data structure 
% clear everything 
[LocFile, LocPath] = uigetfile({'*.mat';'*.*'},'Load in the GREEN data');
if isequal(LocFile,0)
    error('User cancelled the program');
end
load([LocPath LocFile]);

LocFiles = [LocFiles; {[LocPath LocFile]}];

[fidFile, fidPath] = uigetfile({'*.mat';'*.*'},'Load in the RED fid for green data');
if isequal(LocFile,0)
    error('User cancelled the program');
end
load([fidPath fidFile]);

% convert to nm
FidPosR = loc.results.traj.*1e9;
% change coords from relative to absolute
FidPosR_xShift = (min(loc.data.cropCoordsC) + length(loc.data.cropCoordsC)./2).*nmPerPixelR;
FidPosR_yShift = (min(loc.data.cropCoordsR) + length(loc.data.cropCoordsR)./2).*nmPerPixelR;
% rescale relative to starting position of this fiducial
FidPosR(:,1:2) = FidPosR(:,1:2) + [FidPosR_xShift FidPosR_yShift];
FidPosR(:,3) = FidPosR(:,3) - origStartFidPosR;
% avgFidR is the red fiducial for the green data
avgFidR = FidPosR;
% transform red fiducial into green channel
transformedFidRtoG = transformData([avgFidR(:,1),avgFidR(:,2),avgFidR(:,3)],tform);

goodidxnotnan = ~isnan(totalPSFfits(:,zCol));
tempPSFfits = totalPSFfits(goodidxnotnan,:);
clear totalPSFfits
totalPSFfits = tempPSFfits;
clear tempPSFfits

dataSet.frame = totalPSFfits(:,frameCol); 
dataSet.matchConf = totalPSFfits(:,matchConCol); 
dataSet.amp1 = totalPSFfits(:,amp1Col); 
dataSet.amp2 = totalPSFfits(:,amp2Col); 
dataSet.xLoc1 = totalPSFfits(:,xLoc1Col); 
dataSet.yLoc1 = totalPSFfits(:,yLoc1Col); 
dataSet.xLoc2 = totalPSFfits(:,xLoc2Col); 
dataSet.yLoc2 = totalPSFfits(:,yLoc2Col); 
dataSet.sigma1 = totalPSFfits(:,sigma1Col); 
dataSet.sigma2 = totalPSFfits(:,sigma2Col); 
dataSet.bkgndMean = totalPSFfits(:,bkgndMeanCol); 
dataSet.fitError = totalPSFfits(:,fitErrorCol); 
dataSet.goodFitFlag = totalPSFfits(:,goodFitFlagCol); 
dataSet.xCenter = totalPSFfits(:,xCenterCol); 
dataSet.yCenter = totalPSFfits(:,yCenterCol); 
dataSet.angle = totalPSFfits(:,AngleCol); 
dataSet.photons = totalPSFfits(:,PhotonCol); 
dataSet.interlobe = totalPSFfits(:,interlobeCol); 
dataSet.ampRatio = totalPSFfits(:,ampRatioCol); 
dataSet.sigmaRatio = totalPSFfits(:,sigmaRatioCol); 
dataSet.x = totalPSFfits(:,xCol); 
dataSet.y = totalPSFfits(:,yCol); 
dataSet.z = totalPSFfits(:,zCol); 
dataSet.LocFile = LocFile; 
dataSet.LocPath = LocPath; 

dataSet.z = dataSet.z * nSample/nOil; % index mismatch correction

%     dataSet.xFid = totalPSFfits(:,xFidCol); 
%     dataSet.yFid = totalPSFfits(:,yFidCol); 
%     dataSet.zFid = totalPSFfits(:,zFidCol); 
    dataSet.fidTrackX = transformedFidRtoG(:,1); 
    dataSet.fidTrackY = transformedFidRtoG(:,2); 
    dataSet.fidTrackZ = transformedFidRtoG(:,3); 


% transform green data: 
% plot(totalPSFfits(:,xFidCol),totalPSFfits(:,yFidCol),'b.'); 
% 
% transformedData = transformData([totalPSFfits(:,xCol),totalPSFfits(:,yCol),totalPSFfits(:,zCol)],tform); 
% % transformedData = transformData([totalPSFfits(:,xFidCol),totalPSFfits(:,yFidCol),totalPSFfits(:,zFidCol)],tform); 
dataSet.xTrans = NaN; 
dataSet.yTrans = NaN; 
dataSet.zTrans = NaN; 

 % Add the sigmaXYZ for the two colors:
amplitude =  [  361035.867260138,22.2956414971275;...   %   [A1x  A2x]
348907.934759022,28.3183226442783;...   %   [A1y  A2y]
840446.405407229,23.3314294806927];      %   [A1z  A2z]
    
    % Equation 4 of Stallinga and Rieger, ISBI, Barcelona conference proveedings
    
    dataSet.sigmaX = sqrt(amplitude(1,1) .* (1./dataSet.photons) + ...
        amplitude(1,1)*4*amplitude(1,2) .* dataSet.bkgndMean./(dataSet.photons).^2 + ...
        amplitude(1,1) .* (1./dataSet.photons) .* sqrt((2*amplitude(1,2)*(dataSet.bkgndMean./dataSet.photons))./(1+(4*amplitude(1,2)*(dataSet.bkgndMean./dataSet.photons)))));
    
    dataSet.sigmaY = sqrt(amplitude(2,1) .* (1./dataSet.photons) + ...
        amplitude(2,1)*4*amplitude(2,2) .* dataSet.bkgndMean./(dataSet.photons).^2 + ...
        amplitude(2,1) .* (1./dataSet.photons) .* sqrt((2*amplitude(2,2)*(dataSet.bkgndMean./dataSet.photons))./(1+(4*amplitude(2,2)*(dataSet.bkgndMean./dataSet.photons)))));
    
    dataSet.sigmaZ = sqrt(amplitude(3,1) .* (1./dataSet.photons) + ...
        amplitude(3,1)*4*amplitude(3,2) .* dataSet.bkgndMean./(dataSet.photons).^2 + ...
        amplitude(3,1) .* (1./dataSet.photons) .* sqrt((2*amplitude(3,2)*(dataSet.bkgndMean./dataSet.photons))./(1+(4*amplitude(3,2)*(dataSet.bkgndMean./dataSet.photons)))));


% Append the structure to the previous structures in the array
dataSets = [dataSets, dataSet];

clear totalPSFfits avgDevX_denoised avgDevY_denoised avgDevZ_denoised 

%% load in the RED .mat file
% Add it to the total data structure 
% clear everything 

[LocFile, LocPath] = uigetfile({'*.mat';'*.*'},'Load in the RED data');
if isequal(LocFile,0)
    error('User cancelled the program');
end
load([LocPath LocFile]);

LocFiles = [LocFiles; {[LocPath LocFile]}];

[fidFile, fidPath] = uigetfile({'*.mat';'*.*'},'Load in the GREEN fid for RED data');
if isequal(LocFile,0)
    error('User cancelled the program');
end
load([fidPath fidFile]);

% convert to nm
FidPosG = loc.results.traj.*1e9;
% change coords from relative to absolute
FidPosG_xShift = (min(loc.data.cropCoordsC) + length(loc.data.cropCoordsC)./2).*nmPerPixelG;
FidPosG_yShift = (min(loc.data.cropCoordsR) + length(loc.data.cropCoordsR)./2).*nmPerPixelG;
% rescale relative to starting position of this fiducial
FidPosG(:,1:2) = FidPosG(:,1:2) + [FidPosG_xShift FidPosG_yShift];
FidPosG(:,3) = FidPosG(:,3) - origStartFidPosG;
% avgFidR is the red fiducial for the green data
avgFidG = FidPosG;

% avgFidG contains green fiducial for red data

goodidxnotnan = ~isnan(totalPSFfits(:,zCol));
tempPSFfits = totalPSFfits(goodidxnotnan,:);
clear totalPSFfits
totalPSFfits = tempPSFfits;
clear tempPSFfits


dataSet.frame = totalPSFfits(:,frameCol);           % 1
dataSet.matchConf = totalPSFfits(:,matchConCol);    % 2 
dataSet.amp1 = totalPSFfits(:,amp1Col);             % 3
dataSet.amp2 = totalPSFfits(:,amp2Col);             % 4
dataSet.xLoc1 = totalPSFfits(:,xLoc1Col);           % 5
dataSet.yLoc1 = totalPSFfits(:,yLoc1Col);           % 6
dataSet.xLoc2 = totalPSFfits(:,xLoc2Col);           % 7
dataSet.yLoc2 = totalPSFfits(:,yLoc2Col);           % 8
dataSet.sigma1 = totalPSFfits(:,sigma1Col);         % 9
dataSet.sigma2 = totalPSFfits(:,sigma2Col);         % 10 
dataSet.bkgndMean = totalPSFfits(:,bkgndMeanCol);   % 11
dataSet.fitError = totalPSFfits(:,fitErrorCol);     % 12
dataSet.goodFitFlag = totalPSFfits(:,goodFitFlagCol); % 13
dataSet.xCenter = totalPSFfits(:,xCenterCol);         % 14
dataSet.yCenter = totalPSFfits(:,yCenterCol);         % 15
dataSet.angle = totalPSFfits(:,AngleCol);             % 16
dataSet.photons = totalPSFfits(:,PhotonCol);          % 17
dataSet.interlobe = totalPSFfits(:,interlobeCol);     % 18
dataSet.ampRatio = totalPSFfits(:,ampRatioCol);       % 19
dataSet.sigmaRatio = totalPSFfits(:,sigmaRatioCol);   % 20
dataSet.x = totalPSFfits(:,xCol);                     % 21
dataSet.y = totalPSFfits(:,yCol);                     % 22
dataSet.z = totalPSFfits(:,zCol);                     % 23
dataSet.LocFile = LocFile; 
dataSet.LocPath = LocPath; 

%     dataSet.xFid = totalPSFfits(:,xFidCol); 
%     dataSet.yFid = totalPSFfits(:,yFidCol); 
%     dataSet.zFid = totalPSFfits(:,zFidCol); 
    dataSet.fidTrackX = avgFidG(:,1); 
    dataSet.fidTrackY = avgFidG(:,2); 
    dataSet.fidTrackZ = avgFidG(:,3); 



% transform RED data into GREEN channel: 
transformedData = transformData([totalPSFfits(:,xCol),totalPSFfits(:,yCol),totalPSFfits(:,zCol)],tform); 
% transformedData = transformData([totalPSFfits(:,xFidCol),totalPSFfits(:,yFidCol),totalPSFfits(:,zFidCol)],tform); 
dataSet.xTrans = transformedData(:,1); 
dataSet.yTrans = transformedData(:,2); 
dataSet.zTrans = transformedData(:,3); 
dataSet.zTrans = dataSet.zTrans * nSample/nOil; % index mismatch correction

% Equation 4 of Stallinga and Rieger, ISBI, Barcelona conference proveedings
    
    dataSet.sigmaX = sqrt(amplitude(1,1) .* (1./dataSet.photons) + ...
        amplitude(1,1)*4*amplitude(1,2) .* dataSet.bkgndMean./(dataSet.photons).^2 + ...
        amplitude(1,1) .* (1./dataSet.photons) .* sqrt((2*amplitude(1,2)*(dataSet.bkgndMean./dataSet.photons))./(1+(4*amplitude(1,2)*(dataSet.bkgndMean./dataSet.photons)))));
    
    dataSet.sigmaY = sqrt(amplitude(2,1) .* (1./dataSet.photons) + ...
        amplitude(2,1)*4*amplitude(2,2) .* dataSet.bkgndMean./(dataSet.photons).^2 + ...
        amplitude(2,1) .* (1./dataSet.photons) .* sqrt((2*amplitude(2,2)*(dataSet.bkgndMean./dataSet.photons))./(1+(4*amplitude(2,2)*(dataSet.bkgndMean./dataSet.photons)))));
    
    dataSet.sigmaZ = sqrt(amplitude(3,1) .* (1./dataSet.photons) + ...
        amplitude(3,1)*4*amplitude(3,2) .* dataSet.bkgndMean./(dataSet.photons).^2 + ...
        amplitude(3,1) .* (1./dataSet.photons) .* sqrt((2*amplitude(3,2)*(dataSet.bkgndMean./dataSet.photons))./(1+(4*amplitude(3,2)*(dataSet.bkgndMean./dataSet.photons)))));


% Append the structure to the previous structures in the array
dataSets = [dataSets, dataSet];

%% clear all except the dataSets

clearvars -except dataSets LocFiles

% plot no fiducial correction
plot(dataSets(1,2).xTrans,dataSets(1,2).yTrans,'r.')
hold on
plot(dataSets(1,1).x,dataSets(1,1).y,'g.')
axis image 


% APPLY FID corrections

% % green data
% dataSets(1).xFid = dataSets(1).x - dataSets(1).fidTrackX(dataSets(1).frame); 
% dataSets(1).yFid = dataSets(1).y - dataSets(1).fidTrackY(dataSets(1).frame); 
% dataSets(1).zFid = dataSets(1).z - dataSets(1).fidTrackZ(dataSets(1).frame); 
% 
% % red data
% dataSets(2).xFid = dataSets(2).xTrans - dataSets(1).fidTrackX(dataSets(2).frame); 
% dataSets(2).yFid = dataSets(2).yTrans - dataSets(1).fidTrackY(dataSets(2).frame); 
% dataSets(2).zFid = dataSets(2).zTrans - dataSets(1).fidTrackZ(dataSets(2).frame); 

% green data
dataSets(1).xFid = dataSets(1).x - dataSets(1).fidTrackX(dataSets(1).frame); 
dataSets(1).yFid = dataSets(1).y - dataSets(1).fidTrackY(dataSets(1).frame); 
dataSets(1).zFid = dataSets(1).z - dataSets(1).fidTrackZ(dataSets(1).frame); 

% red data
dataSets(2).xFid = dataSets(2).xTrans - dataSets(2).fidTrackX(dataSets(2).frame); 
dataSets(2).yFid = dataSets(2).yTrans - dataSets(2).fidTrackY(dataSets(2).frame); 
dataSets(2).zFid = dataSets(2).zTrans - dataSets(2).fidTrackZ(dataSets(2).frame); 

% plot with fiducial correction
figure; 
plot(dataSets(1,2).xFid,dataSets(1,2).yFid,'r.')
hold on
plot(dataSets(1,1).xFid,dataSets(1,1).yFid,'g.')
axis image 

% plot3(dataSets(1,2).xFid,dataSets(1,2).yFid,dataSets(1,2).zFid,'r.')
% hold on
% plot3(dataSets(1,1).xFid,dataSets(1,1).yFid,dataSets(1,1).zFid,'g.')
% axis image 


%% save the data 
[saveFilePrefix, savePath] = uiputfile({'*.*'},'Enter a prefix for the combined data sets');
if isequal(saveFilePrefix,0)
    error('User cancelled the program');
end
saveFile = [savePath saveFilePrefix];


save([saveFile '.mat'])