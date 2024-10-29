% Copyright (c)2013-2016, The Board of Trustees of The Leland Stanford
% Junior University. All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
% Redistributions of source code must retain the above copyright notice, 
% this list of conditions and the following disclaimer.
% Redistributions in binary form must reproduce the above copyright notice, 
% this list of conditions and the following disclaimer in the documentation 
% and/or other materials provided with the distribution.
% Neither the name of the Leland Stanford Junior University nor the names 
% of its contributors may be used to endorse or promote products derived 
% from this software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS 
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [outputFilePrefix, numBeads,horizontalDHPSF, angledNHA] = ...
    f_calDHPSF(conversionGain,nmPerPixel,boxRadius,channel,sigmaBounds,lobeDistBounds)
% f_calDHPSF is a module in easy_dhpsf that calibrates the z vs. angle
% response of the DH-PSF using one or more discrete fluorescent particles,
% usually fluorescent beads. To generate the data, an objective stepper
% should be used, and the order and size of these steps read out
% from a .dat log file. In addition to the calibration, f_calDHPSF
% generates a series of templates used for template matching.

%% ask user for relevant datafiles
[dataFile, dataPath] = uigetfile({'*.tif';'*.*'},'Open image stack for data processing');
if isequal(dataFile,0)
   error('User cancelled the program');
end
dataFile = [dataPath dataFile];

% hard codes batch processing to be OFF
batchProcess = 0;

% batch processing: turned off unless batchProcess set to 1
if batchProcess == 1
    dlg_title = 'Please Specify Filename Separation for Batch Processing';
    prompt = {  'Filename Common Part:'};
    def = {     dataFile(1:length(dataFile)-4)};
    num_lines = 1;
    inputdialog = inputdlg(prompt,dlg_title,num_lines,def);

    commonRoot = inputdialog{1};

    filelist = dir([fileparts(dataFile) filesep [commonRoot '*.tif']]);  % first characters of file header
    fileNames = {filelist.name}';
else
    fileNames = {dataFile};
end

[darkFile, darkPath] = uigetfile({'*.tif';'*.*'},'Open image stack with dark counts (same parameters as calibration)');
[logFile, logPath] = uigetfile({'*.dat';'*.mat';'*.csv';'*.txt';'*.*'},'Open sequence log file for calibration');
logFile = [logPath logFile];
if isequal(logFile,0)
   error('A sequence log file must be specified for the DHPSF calibration');
end

numFiles = length(fileNames);
dataFileInfo = imfinfo(dataFile);
numFrames = length(dataFileInfo);
imgHeight = dataFileInfo.Height;
imgWidth = dataFileInfo.Width;

% saves in labeled directory if a channel is selected
if channel == '0'
    outputFilePrefix = [dataFile(1:length(dataFile)-4) filesep 'calibration ' ...
        datestr(now,'yyyymmdd HHMM') filesep];

else 
    outputFilePrefix = [dataFile(1:length(dataFile)-4) filesep channel(1) ' calibration ' ...
        datestr(now,'yyyymmdd HHMM') filesep];

end

mkdir(outputFilePrefix);

%% Set Instrument Specific Parameters

% try to grab gain setting from filename assuming form "###G"
% this expression finds numbers followed by 'G' (not followed by letters)
EMGain = regexp(dataFile,'[\d]+(?=(G(?![a-zA-Z])))','match');
if length(EMGain) ~= 1
    EMGain = 300; % i.e. if not in filename as form '...###G'
else
    EMGain = str2num(EMGain{1});
end

dlg_title = 'Set EM Gain, background subtraction type';
prompt = { 'EM Gain (1 if no gain):','Use wavelet subtraction?',...
           'Many NHA holes?','If so, are they at 45 degrees?',...
           'Is DHPSF horizontal when in focus?','Box radius for fit (pix)'}; 
def = {num2str(EMGain),'0','0','1','0',num2str(boxRadius-2)};

num_lines = 1;
inputdialog = inputdlg(prompt,dlg_title,num_lines,def);

EMGain = str2double(inputdialog{1});
if EMGain < 1 || isnan(EMGain)
    warning('EMGain should be a number >= 1. Setting to 1...');
    EMGain = 1;
end

% whether to use wavelet background subtraction, or subtract mean of image
% Off by default to minimize any possible change to shape of
% templates.
useWaveSub = logical(str2num(inputdialog{2}));
fillNHA = logical(str2num(inputdialog{3}));
angledNHA = logical(str2num(inputdialog{4}));
horizontalDHPSF = logical(str2num(inputdialog{5}));
% lower box radius to avoid accidentally fitting bright lobes of two
% adjacent PSFs (problem for NHA with 2.5 um spacing if LA high)
boxRadius = str2num(inputdialog{6});
conversionFactor = conversionGain/EMGain;
ampRatioLimit = 0.5;
sigmaRatioLimit = 0.4;
blurSize = 0.5*160/nmPerPixel;

% Options for lsqnonlin
options = optimset('FunValCheck','on','Diagnostics','off','Jacobian','on','Display','off');
scrsz = get(0,'ScreenSize');

%% Compute dark counts
if ~isequal(darkFile,0)
    darkFile = [darkPath darkFile];
    % Computes average of dark frames for background subtraction
    darkFileInfo = imfinfo(darkFile);
    numDarkFrames = length(darkFileInfo);
    darkAvg = zeros(darkFileInfo(1).Height,darkFileInfo(1).Width);
%     im = tiffread2(darkFile,1,numDarkFrames);
    for frame = 1:numDarkFrames
%         darkAvg = darkAvg + double(im(frame).data);
                        darkAvg = darkAvg + double(imread(darkFile,frame,'Info',darkFileInfo));
    end
    darkAvg = darkAvg/numDarkFrames;
    if ~isequal(size(darkAvg),[imgHeight imgWidth])
        warning('Dark count image and data image stack are not the same size. Resizing dark count image...');
        darkAvg = imresize(darkAvg,[imgHeight imgWidth]);
    end
else
    darkAvg = 0;
end
clear darkFileInfo;

%% Pick region of interest for analysis

% Plots the avg image .tif
dataAvg = zeros(imgHeight,imgWidth);
% im = tiffread2(dataFile,1,min(190,length(dataFileInfo)));
for frame = 1:min(190,length(dataFileInfo))
    dataAvg = dataAvg + double(imread(dataFile,frame,'Info',dataFileInfo));
%     dataAvg = dataAvg + double(im(frame).data));
end
dataAvg = dataAvg/190 - darkAvg;
% clear im

% could change contrast here
hROI=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
imagesc(dataAvg); axis image;colormap hot;
% ROI = imrect(gca,[1 1 128 128]);
if channel == 'g'
    ROI = imrect(gca,[1 1 64 64]);
elseif channel == 'r'
    ROI = imrect(gca,[243 243 64 64]);
else
    ROI = imrect(gca,[1 1 64 64]);
end
title({'Shape box and double-click to choose region of interest for PSF fitting' ...
    mat2str(ROI.getPosition)});
addNewPositionCallback(ROI,@(p) title({'Double-click to choose region of interest for PSF fitting' ...
    ['[xmin ymin width height] = ' mat2str(p,3)]}));
% make sure rectangle stays within image bounds
fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
setPositionConstraintFcn(ROI,fcn);
ROI = round(wait(ROI));
cropHeight = ROI(4);
cropWidth = ROI(3);
close(hROI);

%% Ask user for bead location(s)   
% could change contrast here
hLocs=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
imagesc(dataAvg(ROI(2):ROI(2)+ROI(4)-1, ...
    ROI(1):ROI(1)+ROI(3)-1));axis image;colorbar;colormap hot; 
title('Use LMB to select fiducials. Hit enter to stop or use RMB to select the final fiducial.');
if fillNHA
    title('Select four fiducials that form a large box. Then select fiducials to outline the edge of your "FOV".');
end
hold on;
% User will click on beads in image and then text will be imposed on the
% image corresponding to the bead number.
% moleLocs is a n by 2 array -- x y values for each bead
moleLocs = [];
erase = ones(0,2);
n = 0;
% Loop, collecting bead locations and drawing text to mark them
but = 1;
while but == 1
    [xi,yi,but] = ginput(1);
    if isempty(xi)
        break
    end
    if but == 2;
        erase = [erase; round([xi yi])];
        text(xi,yi,'*','color','red','fontsize',13,'fontweight','bold');
        but = 1;
        continue
    end
    n = n+1;
    text(xi,yi,num2str(n),'color','white','fontsize',13,'fontweight','bold');
    moleLocs(n,:) = round([xi yi]);
end
hold off;

if fillNHA
    moleLocs = f_fill_NHA(moleLocs,nmPerPixel,2500,angledNHA);
    if erase
        for eraseIdx = 1:size(erase,1)
            [~,toErase] = min((moleLocs(:,1)-erase(eraseIdx,1)).^2+(moleLocs(:,2)-erase(eraseIdx,2)).^2);
            moleLocs(toErase,:) = [];
        end
    end
    close(hLocs)
    hLocs=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
        imagesc(dataAvg(ROI(2):ROI(2)+ROI(4)-1, ...
        ROI(1):ROI(1)+ROI(3)-1));axis image;colorbar;colormap hot; 
    for n = 1:length(moleLocs)
        text(moleLocs(n,1),moleLocs(n,2),num2str(n),'color','white','fontsize',13,'fontweight','bold');
    end
end
% moleLocs = round(ginput);
numBeads = size(moleLocs,1);
absLocs = zeros(numBeads,2); % actual positions in nm

saveas(hLocs,[outputFilePrefix 'bead map.png']);
close(hLocs);


%% Initialize data arrays

% [frameNum moleNum amp1 amp2 xMean1 yMean1 xMean2 yMean2 sigma1 sigma2 
%  bkgndMean totalFitError goodFit xCenter yCenter angle numPhotons CFocusPosition,
%  lobeDist ampRatio sigmaRatio]
PSFfits = zeros(numFiles,numFrames*numBeads, 21);
phasemask_position = zeros(numFiles,1);
startTime = tic;

startFrame = zeros(numFiles,200);
endFrame = zeros(numFiles,200);

for n = 1:numFiles

dataFile = fileNames{n};
if numFiles > 1
    phasemask_position(n) = str2num(dataFile(length(commonRoot)+1:length(commonRoot)+1+1));
end
dataFileInfo = imfinfo(dataFile);
numFrames = length(dataFileInfo);


%% Analyze the C-focus scan series

if exist('logFile')
    % either give a list of shutter states and relative z positions, or as
    % a .csv or .dat with one row of frame and one row of stage z positions
    % for PI stage scan, load in .mat file
sifLogData =  importdata(logFile);
if size(sifLogData,2) > 2
    logType = 1; % "siflog" format; assume using averaged steps
    sifLogData(:,4) = sifLogData(:,4).*1000; % conver to nm
    
    useSteps = 1;
    step = 1;               
    buffer = 4;             % advance by two frames to make sure new position is achieved
    i = 2*buffer;

    startFrame(n,1) = i;

    while i <= size(sifLogData,1)-2
        if sifLogData(i+1,1) < 0
            endFrame(n,step) = i;
            step = step + 1;    % current step position
            i = i + buffer;     % advance by buffer frames to make sure new position is achieved
            startFrame(n,step) = i;
        else
            i = i + 1;
        end

    end

    numSteps = step -1;          % the last frames are not part of the scan series
    startFrame_copy = startFrame(n,1:length(endFrame));
    initialPosition = 50.0;

    for i=1:numSteps
        meanStageZ(n,i) = mean(sifLogData(round((startFrame_copy(i)+endFrame(n,i))/2):endFrame(n,i),4));     
        stdStageZ(n,i) = std(sifLogData(round((startFrame_copy(i)+endFrame(n,i))/2:endFrame(n,i)),4));
    end
    clear startFrame_copy

    relativeStepSize = diff(meanStageZ(n,:))';

    for i=1:numSteps
        if i == 1
        meanStageZ(n,i) = initialPosition*1000 ;% in units of nm
        else
        meanStageZ(n,i) = meanStageZ(n,i-1) + relativeStepSize(i-1);% in units of nm
        end
    end
    
    framesToAnalyze = 1:(size(sifLogData,1)-2);
elseif strcmp(logFile(end-2:end), 'csv') 
    logType = 2; % given as frame/abs z. may not be in "steps"
    framesToAnalyze = sifLogData(:,1)';
    stageZPos = sifLogData(:,2); % in nm
    % average positions into steps if staying stationary,
    % allowing for some instrument jitter if stage position is *measured*
    useSteps = min(abs(diff(stageZPos))) < 1;
    meanStageZ(n,:) = stageZPos;
    stdStageZ(n,:) = zeros(size(stageZPos));
    numSteps=max(framesToAnalyze); 
elseif strcmp(logFile(end-2:end), 'mat') 
    logType = 3; % zscan done with PI stage/MATLAB
    load(logFile,'fVec_hat','nFrames_per_step');
    
    fVec_hat = fVec_hat.*-1; % flip sign
    sifLogData = zeros(length(fVec_hat).*nFrames_per_step,2);
    sifLogData(:,1) = 1:length(fVec_hat).*nFrames_per_step;
    
    zVecRep = repmat(fVec_hat, [nFrames_per_step, 1]);
    zVecRep = zVecRep(:);
    sifLogData(:,2) = zVecRep-mean(fVec_hat); %fVec_hat is actual z pos (micron units)
    sifLogData(:,2) = sifLogData(:,2).*1000; % convert to nm
    
    framesToAnalyze = sifLogData(:,1)';
    stageZPos = sifLogData(:,2); % in microns

    startFrame = zeros(1,length(fVec_hat));
    endFrame = startFrame;
    
    for i=1:length(startFrame)
        if i==1
            startFrame(i)=1;
        else
            startFrame(i) = (i-1).*nFrames_per_step + 1;
        end
        
        endFrame(i)=i*nFrames_per_step;
    end
    
    numSteps=length(fVec_hat); 
    for i=1:numSteps
        meanStageZ(n,i) = mean(sifLogData(startFrame(i):endFrame(n,i),2)); %nm
        stdStageZ(n,i) = std(sifLogData(startFrame(i):endFrame(n,i),2));
    end
    
    useSteps = min(abs(diff(stageZPos))) < 1;
    
end

else % if no z log
    error('no z information given!')
end

%% Fit chosen beads throughout entire image stack
h = figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');

for a = 1:size(sifLogData,1)
    
%     im = tiffread2(dataFile,a);
%     dataFrame = double(im.data) - darkAvg;
%     clear im
    dataFrame = double(imread(dataFile,a,'Info',dataFileInfo))-darkAvg;
    data = dataFrame(ROI(2):ROI(2)+ROI(4)-1, ...
    ROI(1):ROI(1)+ROI(3)-1);

    % subtract the background and continue
    if useWaveSub
        bkgndImg = f_waveletBackground(data);
    else
        % subtract median to get rough idea of #photons of fiducial
        bkgndImg = repmat(median(data(:)),size(data));
    end
    data = data - bkgndImg;

    % blur data for more robust peak finding
    dataBlur = imfilter(data, fspecial('gaussian', size(data), blurSize), 'same', 'conv');


    %% do fitting to extract exact locations of DH-PSFs

    % create reconstructed DH-PSF image from fitted data
    reconstructImg = zeros(cropHeight, cropWidth);
    for b=1:numBeads
        rowIdx = (a-1)*numBeads+b;
        
        % create indices to use for fitting
        [xIdx, yIdx] = meshgrid(moleLocs(b,1)-boxRadius:moleLocs(b,1)+boxRadius, ...
            moleLocs(b,2)-boxRadius:moleLocs(b,2)+boxRadius);
        % make sure indices are inside image
        if min(xIdx(:)) < 1
            xIdx = xIdx + (1-min(xIdx(:)));
        end
        if max(xIdx(:)) > cropWidth
            xIdx = xIdx - (max(xIdx(:))-cropWidth);
        end
        if min(yIdx(:)) < 1
            yIdx = yIdx + (1-min(yIdx(:)));
        end
        if max(yIdx(:)) > cropHeight
            yIdx = yIdx - (max(yIdx(:))-cropHeight);
        end
              
        % find two largest peaks in box
        [tempY, tempX] = ind2sub([2*boxRadius+1 2*boxRadius+1], ...
            find(imregionalmax(dataBlur(yIdx(:,1),xIdx(1,:)))));
        tempX = tempX + min(xIdx(:))-1;
        tempY = tempY + min(yIdx(:))-1;
        temp = sortrows([tempX tempY data(sub2ind([cropHeight cropWidth], ...
            tempY,tempX))],-3);
        
        % set initial fitting parameters
        % [amp1 amp2 xMean1 yMean1 xMean2 yMean2 sigma1 sigma2]
        % if two peaks aren't found, then use previous fitting parameters
        if size(temp,1)>=2
            fitParam(3) = temp(1,1);
            fitParam(4) = temp(1,2);
            fitParam(5) = temp(2,1);
            fitParam(6) = temp(2,2);
            fitParam(1) = temp(1,3);
            fitParam(2) = temp(2,3);
            fitParam(7) = mean(sigmaBounds);
            fitParam(8) = mean(sigmaBounds);
        end
        lowerBound = [0 0 min(xIdx(:)) min(yIdx(:)) min(xIdx(:)) min(yIdx(:)) ...
            sigmaBounds(1) sigmaBounds(1)];
        upperBound = [max(max(data(yIdx(:,1),xIdx(1,:)))) ...
            max(max(data(yIdx(:,1),xIdx(1,:)))) ...
            max(xIdx(:)) max(yIdx(:)) max(xIdx(:)) max(yIdx(:)) ...
            sigmaBounds(2) sigmaBounds(2)];

        %% Fit with lsqnonlin
        [fitParam,~,residual,exitflag] = lsqnonlin(@(x) ...
            f_doubleGaussianVector(x,data(yIdx(:,1),xIdx(1,:)),0,xIdx,yIdx),...
            fitParam,lowerBound,upperBound,options);
%         bkgndCnts = bkgndImg(round((fitParam(4)+fitParam(6))/2),round((fitParam(3)+fitParam(5))/2));
        bkgndCnts = 0;
        PSFfits(n,rowIdx,1:13) = [a b fitParam bkgndCnts sum(abs(residual)) exitflag];

        % Calculate midpoint between two Gaussian spots
        % convert from pixels to nm
        PSFfits(n,rowIdx,14) = ((fitParam(3)+fitParam(5))/2)*nmPerPixel;
        PSFfits(n,rowIdx,15) = ((fitParam(4)+fitParam(6))/2)*nmPerPixel;

        % Below is the calculation of the angle of the two lobes.
        % Default is that two vertical lobes is focal plane.
        % Therefore, we want y2>y1 for all angle calculations
        % (so that -90<=angle<=90, and we swap x and y for the atan2
        % calculation. Can also change x/y if default is not correct.
        x1 = fitParam(3);
        x2 = fitParam(5);
        y1 = fitParam(4);
        y2 = fitParam(6);
        % swap if x1>x2
%        if (x1 > x2) 
        if ( ~horizontalDHPSF && (y1 > y2) )  ||  (horizontalDHPSF && x1 > x2)
            tx = x1; ty = y1;
            x1 = x2; y1 = y2;
            x2 = tx; y2 = ty;
            clear tx ty;
        end
        %Finds the angle
        if ~horizontalDHPSF
            PSFfits(n,rowIdx,16) = atan2(-(x2-x1),y2-y1) * 180/pi;
        elseif horizontalDHPSF
            PSFfits(n,rowIdx,16) = atan2(-(y2-y1),x2-x1) * 180/pi;
        end
%        PSFfits(rowIdx,16) = atan2(-(y2-y1),x2-x1) * 180/pi;
        clear x1 x2 y1 y2;

        %Below is a way to count the photons. It integrates the box and
        %subtracts the boxarea*offset from the fit. It is inherently flawed
        %if there happen to be bright pixels inside of the fitting region.
        totalCounts = sum(sum(data(yIdx(:,1),xIdx(1,:))));
        PSFfits(n,rowIdx,17) = totalCounts*conversionFactor;
        
        %The actual position read by the C-focus or PI stage encoder
        if logType == 1 % C focus scan
            PSFfits(n,rowIdx,18) = sifLogData(a,4); 
        elseif logType ==2
            PSFfits(n,rowIdx,18) = meanStageZ(a); 
        elseif logType == 3 % PI stage scan
            PSFfits(n,rowIdx,18) = sifLogData(a,2);
        end
        %The interlobe distance
        lobeDist = sqrt((fitParam(3)-fitParam(5)).^2 + ...
                        (fitParam(4)-fitParam(6)).^2);
        PSFfits(n,rowIdx,19) = lobeDist;
        
        %Amplitude Ratio
        ampRatio = abs(fitParam(1) - fitParam(2))/sum(fitParam(1:2));
        PSFfits(n,rowIdx,20) = ampRatio;
        
        % Gaussian width Ratio
        sigmaRatio = abs(fitParam(7) - fitParam(8))/sum(fitParam(7:8));
        PSFfits(n,rowIdx,21) = sigmaRatio;

        %% Now evaluate the fits
        % Conditions for fits (play with these):
        % (1) Amplitude of both lobes > 0
        % (2) All locations x1,y1, x2,y2 lie inside area of small box
        % (3) All sigmas need to be > sigmaBound(1) and < sigmaBound(2)
        % (4) Distance between lobes needs to be > lobeDist(1) pixels and < lobeDist(2) pixels
        % (5) Make sure amplitudes are within 100% of one another
        % (6) Make sure totalFitError/(total number of photons) < 1.05 (not
        %     enabled at the present time)
        if exitflag > 0
            if fitParam(1)<0 || fitParam(2)<0 
                PSFfits(n,rowIdx,13) = -1001;   %% this flag will never show up because it is equivalent to extiflag -2 in lsqnonlin
            end
            if fitParam(3)<min(xIdx(:)) || fitParam(5)<min(xIdx(:)) ...
               || fitParam(3)>max(xIdx(:)) || fitParam(5)>max(xIdx(:)) ...
               || fitParam(4)<min(yIdx(:)) || fitParam(6)<min(yIdx(:)) ...
               || fitParam(4)>max(yIdx(:)) || fitParam(6)>max(yIdx(:))
                PSFfits(n,rowIdx,13) = -1002;   %% this flag will never show up because it is equivalent to extiflag -2 in lsqnonlin
            end
            if fitParam(7)<=sigmaBounds(1) || fitParam(8)<=sigmaBounds(1) ...
                || fitParam(7)>=sigmaBounds(2) || fitParam(8)>=sigmaBounds(2)
                PSFfits(n,rowIdx,13) = -1003;   
            end
            if sigmaRatio > sigmaRatioLimit;
                PSFfits(n,rowIdx,13) = -1004;
            end
            if lobeDist < lobeDistBounds(1) || lobeDist > lobeDistBounds(2)        
                PSFfits(n,rowIdx,13) = -1005;
            end
            if ampRatio > ampRatioLimit;
                PSFfits(n,rowIdx,13) = -1006;
            end
%             if PSFfits(n,rowIdx,12)*conversionFactor/PSFfits(n,rowIdx,17) > 1 || ...
%                PSFfits(n,rowIdx,12)*conversionFactor/PSFfits(n,rowIdx,17) < 0
%                 PSFfits(n,rowIdx,13) = -1007;
%             end
            % does the fit jump a large amount (>400 nm) from previous fit?
            if norm(squeeze(PSFfits(n,rowIdx,14:15))'-moleLocs(b,:)*nmPerPixel) > 400;
                PSFfits(n,rowIdx,13) = -1008;
            end
        end
        
        % if fit was successful, use the computed center location as center
        % of box for next iteration
        if PSFfits(n,rowIdx,13) > 0
            moleLocs(b,:) = round(PSFfits(n,rowIdx,14:15)/nmPerPixel);            
        end
        
        % plot image reconstruction so that fits can be checked
        [xIdx, yIdx] = meshgrid(1:cropWidth,1:cropHeight);
        reconstructImg = reconstructImg + ...
            fitParam(1).*exp( -((xIdx-fitParam(3)).^2+(yIdx-fitParam(4)).^2.) / (2.*fitParam(7).^2)) ...
            +fitParam(2).*exp( -((xIdx-fitParam(5)).^2+(yIdx-fitParam(6)).^2.) / (2.*fitParam(8).^2));
    end

%%  plot results of template matching and fitting
    set(0,'CurrentFigure',h);
    
    % could change contrast here
    data = data+bkgndImg;
    subplot('Position',[0.025 0.025 .9/2 .95]);
    imagesc(data);axis image;colormap hot;
    title(['Frame ' num2str(a) ': raw data - dark counts']);

    subplot('Position',[0.525 0.025 .9/2 .95]);
    imagesc(reconstructImg+bkgndImg,[min(data(:)),max(data(:))]);
    axis image;
    title('Image reconstructed from fitted matches');

    drawnow;
end

elapsedTime = toc(startTime);
clear data residual dataAvg reconstructImg xIdx yIdx temp;
close(h);
fps = numFrames/elapsedTime
beadsPerSec = numFrames*numBeads/elapsedTime

for b = 1:numBeads
    absLocs(b,1:2) = mean(squeeze(PSFfits(1,PSFfits(1,:,2)==b,14:15)),1) + [ROI(1)-1,ROI(2)-1] * nmPerPixel; 
end

% save fit info to MATLAB mat file
save([outputFilePrefix 'raw fits.mat']);


end

%% initialize calibration variables

% From the old version, used when only using one file and one bead
%numSteps = numFrames/numAcqPerStep;
% meanAngles = zeros(1,numSteps);
% stddevAngles = zeros(1,numSteps);
% meanPhotons = zeros(1,numSteps);
% stddevPhotons = zeros(1,numSteps);
% meanX = zeros(1,numSteps);
% stdX = zeros(1,numSteps);
% meanY = zeros(1,numSteps);
% stdY = zeros(1,numSteps);
% numGoodFrames = zeros(1,numSteps);

% load('raw fits.mat')

meanAngles = zeros(numFiles,numBeads,numSteps);
stddevAngles = zeros(numFiles,numBeads,numSteps);
meanPhotons = zeros(numFiles,numBeads,numSteps);
stddevPhotons = zeros(numFiles,numBeads,numSteps);
meanX = zeros(numFiles,numBeads,numSteps);
stdX = zeros(numFiles,numBeads,numSteps);
meanY = zeros(numFiles,numBeads,numSteps);
stdY = zeros(numFiles,numBeads,numSteps);
numGoodFrames = zeros(numFiles,numBeads,numSteps);

meanInterlobeDistance = zeros(numFiles,numBeads,numSteps);
stdInterlobeDistance = zeros(numFiles,numBeads,numSteps);
meanAmp1 = zeros(numFiles,numBeads,numSteps);
stdAmp1 = zeros(numFiles,numBeads,numSteps);
meanAmp2 = zeros(numFiles,numBeads,numSteps);
stdAmp2 = zeros(numFiles,numBeads,numSteps);
meanAmpRatio = zeros(numFiles,numBeads,numSteps);
stdAmpRatio = zeros(numFiles,numBeads,numSteps);
meanSigma1 = zeros(numFiles,numBeads,numSteps);
stdSigma1 = zeros(numFiles,numBeads,numSteps);
meanSigma2 = zeros(numFiles,numBeads,numSteps);
stdSigma2 = zeros(numFiles,numBeads,numSteps);
meanSigmaRatio = zeros(numFiles,numBeads,numSteps);
stdSigmaRatio = zeros(numFiles,numBeads,numSteps);
goodFit_f = zeros(numFiles,numBeads,numSteps);
goodFit_b = zeros(numFiles,numBeads,numSteps);
zAngleZero = zeros(numFiles,numBeads);
z = zeros(numFiles,numBeads,numSteps);
%     figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');

textHeader = {'Start Frame' 'End Frame' 'Mean Angle (deg)' ...
    'x deviation (nm)' 'y deviation (nm)' 'z position (nm)' 'Mean number of Photons' ...
    'Angle Std Dev (deg)' 'x std dev (nm)' 'y std dev (nm)' 'z Std Dev (nm)' ...
    'Photons Std Dev' 'Num good frames'};

% compute calibration parameters
for n = 1:numFiles
    for bead = 1:numBeads
        
    % extract fitting parameters for this bead
    beadFitParam(:,:) = PSFfits(n, PSFfits(n,:,2) == bead, :);

    for step = 1:numSteps
%     for step = 8:18
        % extract bead fitting parameters for this step
        
%             if numFrames ~= numSteps 
%                 stepFitParam = beadFitParam(startFrame(n,step):endFrame(n,step),:);
%             else
%                 stepFitParam = beadFitParam(step,:);
%             end

            stepFitParam = beadFitParam(startFrame(n,step):endFrame(n,step),:);

            
        x = stepFitParam(:,14);
        y = stepFitParam(:,15);
        angles = stepFitParam(:,16);

        % correct noise in angle measurement if oscillating between +90
        % and -90 degrees
        if (~isempty(angles(angles > 80)) && ~isempty(angles(angles < -80)))
%            fprintf(2,[mat2str(angles) '\n']);
            angles(angles < 0) = angles(angles < 0) + 180;
%            fprintf(2,[mat2str(angles) '\n']);
        end

        amp1 = stepFitParam(:,3);
        amp2 = stepFitParam(:,4);
        sigma1 = stepFitParam(:,9);
        sigma2 = stepFitParam(:,10);
        numPhotons = stepFitParam(:,17);
        lobeDist = stepFitParam(:,19);
        ampRatio = stepFitParam(:,20);
        sigmaRatio = stepFitParam(:,21);

        goodFit = stepFitParam(:,13)>0;

        meanX(n, bead,step) = mean(x(goodFit));
        stdX(n,bead,step) = std(x(goodFit));
        meanY(n,bead,step) = mean(y(goodFit));
        stdY(n,bead,step) = std(y(goodFit));
        meanAngles(n,bead,step) = mean(angles(goodFit));
        stddevAngles(n,bead,step) = std(angles(goodFit));
        meanPhotons(n,bead,step) = mean(numPhotons(goodFit));
        stddevPhotons(n,bead,step) = std(numPhotons(goodFit));
        numGoodFrames(n,bead,step) = length(angles(goodFit));
        
        meanInterlobeDistance(n,bead,step) = mean(lobeDist(goodFit));
        stdInterlobeDistance(n,bead,step) = std(lobeDist(goodFit));
        meanAmp1(n,bead,step) = mean(amp1(goodFit));
        stdAmp1(n,bead,step) = std(amp1(goodFit));
        meanAmp2(n,bead,step) = mean(amp2(goodFit));
        stdAmp2(n,bead,step) = std(amp2(goodFit)); 
        meanSigma1(n,bead,step) = mean(sigma1(goodFit));
        stdSigma1(n,bead,step) = std(sigma1(goodFit));
        meanSigma2(n,bead,step) = mean(sigma2(goodFit));
        stdSigma2(n,bead,step) = std(sigma2(goodFit));
        meanAmpRatio(n,bead,step) = mean(ampRatio(goodFit));
        stdAmpRatio(n,bead,step) = std(ampRatio(goodFit));       
        meanSigmaRatio(n,bead,step) = mean(sigmaRatio(goodFit));
        stdSigmaRatio(n,bead,step) = std(sigmaRatio(goodFit));       

    end

    
    %unwrap angles if we have greater than 180 degree range; need to convert
    %180 degrees to 2pi because unwrap function works on intervals of 2pi, not
    %pi
    meanAngles(n,bead,:) = unwrap(meanAngles(n,bead,:)*2*pi/180)*180/(2*pi);
    if isempty(meanAngles(n,bead,meanAngles(n,bead,:)>10))
        meanAngles(n,bead,:) = meanAngles(n,bead,:) + 180;
    end
    if isempty(meanAngles(n,bead,meanAngles(n,bead,:)<-10))
        meanAngles(n,bead,:) = meanAngles(n,bead,:) - 180;
    end


    % this computes the differences in angle between adjacent steps
    % these differences are sampled at half points between steps
    angleSlope = squeeze(squeeze(diff(meanAngles(n,bead,:))));
    % average adjacent differences to resample differences at integer points,
    % convert to units of angle/nm
    indices = 1:length(angleSlope)-1;
    stepSize =  squeeze(diff(meanStageZ(n,:)));           %in units of nm
    anglePerNM = (angleSlope(indices+1)+angleSlope(indices))./(2*stepSize(indices)');  %check if correct
    stddevNM = squeeze(squeeze(stddevAngles(n,bead, 2:length(stddevAngles(n,bead,:))-1)))./anglePerNM*sign(median(anglePerNM,'omitnan'));

    % test to make sure angle vs z curve is monotonic -- remove any
    % outlying points
%     slope = sign(meanAngles(n,bead,length(meanAngles(n,bead,:)))-meanAngles(n,bead,1));
    slope = sign(meanAngles(n,bead,round(numSteps*1/4)+3)-meanAngles(n,bead,round(numSteps*1/4)-3));
    
    goodFit = squeeze(squeeze(~isnan(meanAngles(n,bead,:))))';
    
    stepSlope = false(1,length(angleSlope));
    indices_forward = find(sign(angleSlope)==1);
    stepSlope(indices_forward)=true;
    goodFit_forward = logical([0 stepSlope]) & goodFit;
    stepSlope = false(1,length(angleSlope));
    indices_backward = find(sign(angleSlope)==-1);
    stepSlope(indices_backward(2:length(indices_backward)))=true;
    goodFit_backward = logical([0 stepSlope]) & goodFit;
    
    % if calibration movie only contains one backward scan,
    % then replicate it as if it were moving forward
    % Try to use forward frames unless only a few are fit, or many more
    % backward are fit.
    if sum(goodFit_backward) > sum(goodFit_forward)*2
        gfTemp = goodFit_forward;
        goodFit_forward = goodFit_backward;
        goodFit_backward = gfTemp;
        clear gfTemp
    end

%     goodFit = logical([0 ones(1,13) zeros(1,11)]) & ...
%         squeeze(squeeze(~isnan(meanAngles(n,bead,:))))' & ...
%         [true squeeze(sign(meanAngles(n,bead,2:size(meanAngles,3))-...
%         meanAngles(n,bead,1:size(meanAngles,3)-1)))'==slope];

%     goodFit = goodFit(2:length(goodFit)-1);
%     z = 0:stepSize:stepSize*(numSteps-1);
    if logType == 1
        z0 = meanStageZ(n,:) - meanStageZ(n,1); % assuming begins in focus
    elseif logType ==2
        z0 = meanStageZ(n,:);
    elseif logType == 3
        z0 = meanStageZ(n,:); %%%%
    end
    z1 = z0(2:length(z0)-1);
    % flip definition of z if angle vs z slope is positive (should be
    % negative)
    if slope>0
        z0 = -z0;
        z1 = -z1;
    end

    % if there aren't enough points for a good curve, skip this fiduciary
    if sum(goodFit_forward) < 5
        disp(['Fiducial number ' num2str(bead) ' had only ' num2str(sum(goodFit_forward)) ' usable steps and is excluded']);
        continue;
    end
    if logType == 1 || logType == 3% as usual, define relative space
    % compute xyz vector so that z = 0 corresponds to angle = 0
    xAngleZero = interp1(squeeze(meanAngles(n,bead,goodFit_forward)),...
        squeeze(meanX(n,bead,goodFit_forward)),0,'spline');
    yAngleZero = interp1(squeeze(meanAngles(n,bead,goodFit_forward)),...
        squeeze(meanY(n,bead,goodFit_forward)),0,'spline');
    zAngleZero(n,bead) = interp1(squeeze(meanAngles(n,bead,goodFit_forward)),...
        z0(goodFit_forward),0,'spline');
    elseif logType ==2 % defined for contest
        xAngleZero = meanX(n,bead,76);
        yAngleZero = meanY(n,bead,76);
        zAngleZero(n,bead) = 0;
        goodFit_forward = goodFit;
        goodFit_backward = zeros(size(goodFit));
    end
    meanX(n,bead,goodFit) = meanX(n,bead,goodFit) - xAngleZero;
    meanY(n,bead,goodFit) = meanY(n,bead,goodFit) - yAngleZero;
    z0(goodFit) = z0(goodFit)-zAngleZero(n,bead);
    z1 = z1-zAngleZero(n,bead);

    
    %% plot calibration parameters
    if logical(sum(goodFit_backward))
        scanLegend = {'forward scan','backward scan'};
        precLegend = {'x forward','x backward','y forward','y backward'};
    else
        scanLegend = {'calibration scan'};
        precLegend = {'x precision','y precision'};
    end
    h1 = figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
    subplot(2,3,1:2);
    if logType == 1 || logType == 3 % any(stddevAngles(n,bead,goodFit_forward))
        errorbar(z0(goodFit_forward),squeeze(meanAngles(n,bead,goodFit_forward)),squeeze(stddevAngles(n,bead,goodFit_forward)),'-');
    elseif logType == 2
        plot(z0(goodFit_forward),squeeze(meanAngles(n,bead,goodFit_forward)),'.-');
    end
    if logical(sum(goodFit_backward))
        hold on
        errorbar(z0(goodFit_backward),squeeze(meanAngles(n,bead,goodFit_backward)),squeeze(stddevAngles(n,bead,goodFit_backward)),':');
        hold off
    end
    axis tight;
    title({dataFile ['Fiduciary ' num2str(bead)]},'interpreter','none');
    legend(scanLegend);
    xlabel('z position (nm)');
    ylabel('Angle (deg)');
    
    subplot(2,3,3);
    plot(z0(goodFit_forward),squeeze(squeeze(meanX(n,bead,goodFit_forward))),'b');
    hold on;
    plot(z0(goodFit_forward),squeeze(squeeze(meanY(n,bead,goodFit_forward))),'r');
    if logical(sum(goodFit_backward))
        plot(z0(goodFit_backward),squeeze(squeeze(meanX(n,bead,goodFit_backward))),':b');
        plot(z0(goodFit_backward),squeeze(squeeze(meanY(n,bead,goodFit_backward))),':r');
    end
    axis tight;
    %     legend('x forward','x backward','y forward','y backward');
    xlabel('z position (nm)');
    ylabel('xy position (nm)');
    hold off;
    
    subplot(2,3,4);
    plot(z0(goodFit_forward),squeeze(stdX(n,bead,goodFit_forward)),'b');
    hold on;
    plot(z0(goodFit_forward),squeeze(stdY(n,bead,goodFit_forward)),'r');
    if logical(sum(goodFit_backward))
        plot(z0(goodFit_backward),squeeze(stdX(n,bead,goodFit_backward)),':b');
        plot(z0(goodFit_backward),squeeze(stdY(n,bead,goodFit_backward)),':r');
    end
    axis tight;
    ylim([0 40]);
    legend(precLegend);
    xlabel('z position (nm)');
    ylabel('localization precision (nm)');
    hold off;
    
    subplot(2,3,5);
    plot(z1(goodFit_forward(2:length(goodFit_forward)-1)),stddevNM(goodFit_forward(2:length(goodFit_forward)-1)),'b');
    if logical(sum(goodFit_backward))
        hold on
        plot(z1(goodFit_backward(2:length(goodFit_backward)-1)),stddevNM(goodFit_backward(2:length(goodFit_backward)-1)),':b');
        hold off
    end
    axis tight;
    ylim([0 40]);
    xlabel('z position (nm)');
    ylabel('z localization precision (nm)');
    
    subplot(2,3,6);
    errorbar(z0(goodFit_forward),squeeze(meanPhotons(n,bead,goodFit_forward)),squeeze(stddevPhotons(n,bead,goodFit_forward)),'b');
    if logical(sum(goodFit_backward))
        hold on
        errorbar(z0(goodFit_backward),squeeze(meanPhotons(n,bead,goodFit_backward)),squeeze(stddevPhotons(n,bead,goodFit_backward)),':b'); 
        hold off
    end
    axis tight;
    xlabel('z position (nm)');
    ylabel('number of photons');
    drawnow
    print(h1,'-dpng',[outputFilePrefix 'bead ' num2str(bead) ' stats_1.png']);
    % close this window if there are a lot of beads that were fitted
    if numBeads > 10
        close(h1);
    end
    
    
    %% plot additional pertaining to the DH-PSF
    
    h2=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
    subplot(2,3,1);
    errorbar(z0(goodFit_forward),squeeze(meanAngles(n,bead,goodFit_forward)),squeeze(stddevAngles(n,bead,goodFit_forward)),'-');
    hold on
%     errorbar(z(goodFit_backward),meanAngles(n,bead,goodFit_backward),stddevAngles(n,bead,goodFit_backward),':');
    axis tight;
    title({dataFile ['Fiduciary ' num2str(bead)]},'interpreter','none');
    legend('forward scan'); %,'backward scan');
    xlabel('z position (nm)');
    ylabel('Angle (deg)');
    hold off

    subplot(2,3,4);
    errorbar(z0(goodFit_forward),squeeze(meanInterlobeDistance(n,bead,goodFit_forward)),squeeze(stdInterlobeDistance(n,bead,goodFit_forward)),'-');
    axis tight;
    xlabel('z position (nm)');
    ylabel('Interlobe Distance (pix)');
    hold off;
    
    subplot(2,3,2);
    errorbar(z0(goodFit_forward),squeeze(meanAmp1(n,bead,goodFit_forward)),squeeze(stdAmp1(n,bead,goodFit_forward)),'-');
    hold on;
    errorbar(z0(goodFit_forward),squeeze(meanAmp2(n,bead,goodFit_forward)),squeeze(stdAmp2(n,bead,goodFit_forward)),'-r');
    axis tight;
    legend('Lobe 1','Lobe 2');
    xlabel('z position (nm)');
    ylabel('Gaussian Amplitudes (counts)');
    hold off;

    subplot(2,3,5);
    errorbar(z0(goodFit_forward),squeeze(meanAmpRatio(n,bead,goodFit_forward)),squeeze(stdAmpRatio(n,bead,goodFit_forward)),'-');
    hold on;
    axis tight;
    xlabel('z position (nm)');
    ylabel('Gaussian Amplitude Ratio');
    hold off;
    
    subplot(2,3,3);
    errorbar(z0(goodFit_forward),squeeze(meanSigma1(n,bead,goodFit_forward)),squeeze(stdSigma1(n,bead,goodFit_forward)),'-');
    hold on;
    errorbar(z0(goodFit_forward),squeeze(meanSigma2(n,bead,goodFit_forward)),squeeze(stdSigma2(n,bead,goodFit_forward)),'-r');
    axis tight;
    legend('Lobe 1','Lobe 2');
    xlabel('z position (nm)');
    ylabel('Gaussian Sigma (pixel)');
    hold off;

    subplot(2,3,6);
    errorbar(z0(goodFit_forward),squeeze(meanSigmaRatio(n,bead,goodFit_forward)),squeeze(stdSigmaRatio(n,bead,goodFit_forward)),'-');
    hold on;
    axis tight;
    xlabel('z position (nm)');
    ylabel('Gaussian Sigma Ratio');
    hold off;
    
    print(h2,'-dpng',[outputFilePrefix 'bead ' num2str(bead) ' stats_2.png']);
    close(h2);
  
    %% move 'local' calibration parameters to 'global' array
    goodFit_f(n,bead,:) = goodFit_forward;
    goodFit_b(n,bead,:) = goodFit_backward;
    z(n,bead,:) = z0;
    
%     bead
    end
% n
end

% average beads together to be slightly more accurate
if isequal(logType,2)
    meanX = mean(meanX,2);
    meanY = mean(meanY,2);
    meanAngles = mean(meanAngles,2);
    disp('Averaging x,y, and z calibration between beads!');
end

%% write calibration parameters to a file
save([outputFilePrefix 'calibration.mat'], ...
        'meanAngles', 'meanX', 'meanY', 'z','xAngleZero','yAngleZero',...
        'zAngleZero', 'goodFit_f', 'absLocs',...
        'goodFit_b', 'meanPhotons', 'stddevPhotons','stdX','stdY',...
        'stddevAngles','meanInterlobeDistance','stdInterlobeDistance',...
        'meanAmpRatio','stdAmpRatio','meanAmp1','meanAmp2','horizontalDHPSF');
    
%% Generate template Stack of a chosen bead 

templateSize = 20;      % 2*round(10*160/nmPerPixel);  % 26;

% make sure templateSize is an even number of pixels
% if mod(templateSize,2)==1
%     templateSize = templateSize+1;
% end



for bead = 1:numBeads
    goodFit_forward = logical(squeeze(goodFit_f(1,bead,:)));
    numTemplates=length(startFrame(goodFit_forward));
   
    if numFrames ~= numSteps 
        forwardStartFrames = startFrame(goodFit_forward);
        forwardEndFrames = endFrame(goodFit_forward);
    else
       	forwardStartFrames=framesToAnalyze(goodFit_forward);
        forwardEndFrames = framesToAnalyze(goodFit_forward);
    end

    template = zeros(numTemplates,templateSize,templateSize);
    for a=1:numTemplates
        for b = round((forwardStartFrames(a)+forwardEndFrames(a))/2):forwardEndFrames(a)

            % load data frames
            
%             im = tiffread2(dataFile,b);
%             dataFrame = double(im.data) - darkAvg;
%             clear im
            dataFrame = double(imread(dataFile,b,'Info',dataFileInfo))-darkAvg;
            data = dataFrame(ROI(2):ROI(2)+ROI(4)-1, ...
            ROI(1):ROI(1)+ROI(3)-1);

            % subtract the background and continue
            if useWaveSub
                bkgndImg = f_waveletBackground(data);
                data = data - bkgndImg;    
            end
            
            good = PSFfits(n,:,1)==b & PSFfits(n,:,2)==bead;
            x_Pos = round(PSFfits(n,good,14)/nmPerPixel);
            y_Pos = round(PSFfits(n,good,15)/nmPerPixel);

            % check to make sure x_Pos,y_Pos are valid values
            if y_Pos-floor(templateSize/2) < 1
                y_Pos = floor(templateSize/2)+1;
            end
            if y_Pos+floor(templateSize/2) > size(data,1)
                y_Pos = size(data,1)-floor(templateSize/2);
            end
            if x_Pos-floor(templateSize/2) < 1
                x_Pos = floor(templateSize/2)+1;
            end
            if x_Pos+floor(templateSize/2) > size(data,2)
                x_Pos = size(data,2)-floor(templateSize/2);
            end
            DHPSF_Image = data(y_Pos-floor(templateSize/2):y_Pos+floor(templateSize/2),...
                x_Pos-floor(templateSize/2):x_Pos+floor(templateSize/2));
            DHPSF_Image = DHPSF_Image(1:templateSize,1:templateSize);   % resize to the correct dimension if necessary
            template(a,:,:) = squeeze(template(a,:,:)) + DHPSF_Image;

        end
    %     imagesc(squeeze(template(a,:,:)))
    %     axis square
    %     drawnow
    %     pause(1)


    end

    save([outputFilePrefix 'bead ' num2str(bead) ' templates.mat'], 'template')

end
end