% Copyright (c)2013, The Board of Trustees of The Leland Stanford Junior
% University. All rights reserved.
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

% other approach to passing variables to this function
% function f_processFits(totalPSFfits,numFrames,ROI,conversionFactor,...
%     sigmaBounds,lobeDistBounds,ampRatioLimit,sigmaRatioLimit,nmPerPixel)
% current approach: load all relevant variables from f_fitSMs output
% directly
function f_processFits(catPSFfits,numFrames,fitFilePrefix, fidTrackX, fidTrackY,...
    fidTrackZ, nmPerPixel,spatialCorr,useCurrent,currFidIdx,fidsToUse,fidChannel,imgShiftX,imgShiftY)
useTimeColors = 0;
plotAsTracks = 0;
showTemplates = 0;
numPhotonRange = [300 10000000];
xyPrecRange = [0 150];
zPrecRange = [0 200];

numFramesAll = sum(numFrames);
load([fitFilePrefix{1} 'molecule fits.mat']);

threshVals = peakThreshold(1,:);

% Parameters

% dlg_title = 'Please Input Parameters';
% prompt = {  'Use Fiducial correction?'};    %,...
% %     'Size of Points in reconstruction',...
% %     'White Light Shift X (in nm)',...
% %     'White Light Shift Y (in nm)',...
% %    'Laser Power at objective (in mW)'...
% %     };
% def = {    '1'};    %, ... 
% %     '30', ...
% %     '0', ...
% %     '0', ...
% %     'NaN' ...
% %     };
% num_lines = 1;
% inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
% 
% useFidCorrections = str2double(inputdialog{1});
scatterSize = 30; %str2double(inputdialog{2});
wlShiftX = 0; %str2double(inputdialog{3});
wlShiftY = 0; %str2double(inputdialog{4});
% powerAtObjective = str2double(inputdialog{5})/1000;
    if sum(isnan(catPSFfits(:,30))) == size(catPSFfits,1)
        warning(['No fiducial correction applied!']);
    end

%% define plotting parameters
% useFidCorrections = logical(useFidCorrections);

% nmPerPixel = 125.78;    % was 160 for 8b back
% scaleBarLength = 1000;  % nm
% pixelSize = 2;          % size of pixels in reconstructed image in nm
% border = 500;           % plot extra region around the cells (size of extra region in nm)
%wlShiftX = 0;          % shift white light image in x direction (in nm)
% (positive = move right)
%wlShiftY = 0;         	% shift white light image in y direction (in nm)
% (positive = move up)
% lambda = 615;           % nm, was 527
% NA = 1.4;               % numerical aperture
nSample = 1.33;         % index of refraction of sample
nOil = 1.518;           % index of immersion oil

scrsz = get(0,'ScreenSize');
% c_map = hot(256);

%% open datafiles

% [locFile, locPath] = uigetfile({'*.mat';'*.*'},'Open data file with PSF localizations');
% if isequal(locFile,0)
%     error('User cancelled the program');
% end
[whiteLightFile, whiteLightPath] = uigetfile({'*.tif';'*.*'},'Open image stack with white light image');
whiteLightFile = [whiteLightPath whiteLightFile];
% load([locPath locFile]);



ROI_initial = ROI + [imgShiftX imgShiftY 0 0]; % assuming WL ROI is full chip

% filePrefix = [dataPath dataFile(1:length(dataFile)-4) ' ' datestr(now,'yyyymmdd HHMM')];


%% Evaluate the laser intensity from the Background fits

% [laser_x_nm, laser_y_nm ,sigma_x_nm, sigma_y_nm, theta, peakIntensity, waist]...
%     = EstimateGaussianLaserProfile...
%     (bkgndImg_avg, FOWmask, nmPerPixel, powerAtObjective);


% laserAmpTreshold = 20;
% 
% laserX =  (mean(bkgndFits(bkgndFits(:,2)>laserAmpTreshold,3))+ROI_initial(1)) * nmPerPixel;            % laser center position
% laserY =  (mean(bkgndFits(bkgndFits(:,2)>laserAmpTreshold,4))+ROI_initial(2)) * nmPerPixel;
% laserWidthX = 4 * mean(bkgndFits(bkgndFits(:,2)>laserAmpTreshold,5)) * nmPerPixel;    % Gaussian Beam width, the factor of 4 is needed to convert
% laserWidthY = 4 * mean(bkgndFits(bkgndFits(:,2)>laserAmpTreshold,6)) * nmPerPixel;    % for the definition of the fitting function to Gaussian intensity function
% laserRot = mean(bkgndFits(bkgndFits(:,2)>laserAmpTreshold,7));
% 
% % Assuming a radially symmetric Gaussian intensity distribution
% peakIntensity = 4*(powerAtObjective)/(2*pi*((laserWidthX+laserWidthY)/2)^2)...
%     * (10^7)^2  % in units of Watts/cm^2

%% Plot reconstructions

% totalPSFfits_original = totalPSFfits;       % copy the original data to avoid corruption
pass = 1;
anotherpass = true;

% if ~exist('numFrames', 'var')
%     numFrames = frames(length(frames));
% end

zRange = [-2000 2000];
%if exist('numFramesAll');%TODO
frameRange = [1 sum(numFramesAll)];
fitErrorRange = [0 4];

while anotherpass == true


    %% Plot the filter parameters
    
    if pass == 1
        
        fitErrorCol = 16;
        goodFitFlagCol = 17;
        numPhotonCol = 21;
        lobeDistCol = 22;
        ampRatioCol = 23;
        sigmaRatioCol = 24;
        templateNumCol = 5;
        templateStrCol = 6;
        
        %totalPSFfits=[frame#, loc# w/in frame, xLocation yLocation matchingTemplateNumber matchConfidence (6),
        %amp1 amp2 xMean1 yMean1 xMean2 yMean2 sigma1 sigma2 bkgndMean(15)
        %totalFitError goodFit xCenter yCenter angle numPhotons
        %interlobeDistance, amplitude ratio, sigma ratio(24)
        %corrected X, corrected Y, corrected Z (27)]
        
        initGoodFits = catPSFfits(:,goodFitFlagCol) > 0;
        initGoodFits = ~isnan(catPSFfits(:,25)); %%%%
        
        
        figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
        subplot(2,2,1)
        [n,xout] = hist(catPSFfits(:,lobeDistCol), 4:0.1:12);
        bar(xout,n)
        hold on
        [n,xout] = hist(catPSFfits(initGoodFits,lobeDistCol), 4:0.1:12);
        bar(xout,n, 'green')
%         [n,xout] = hist(catPSFfits(~initGoodFits,lobeDistCol), 4:0.1:12);
%         bar(xout,n, 'red')
        hold off
        xlabel('pixel'); ylabel('Frequency');
        title('Lobe Distance');
        legend('unfiltered','initially included')
        xlim([4 14]);
        
        subplot(2,2,2)
        unfilteredFitError = catPSFfits(:,fitErrorCol)*conversionFactor./catPSFfits(:,numPhotonCol);        
        fitError = catPSFfits(initGoodFits,fitErrorCol)*conversionFactor./catPSFfits(initGoodFits,numPhotonCol);
        badFitError = catPSFfits(~initGoodFits,fitErrorCol)*conversionFactor./catPSFfits(~initGoodFits,numPhotonCol);
        [n,xout] = hist(unfilteredFitError(unfilteredFitError > 0 & unfilteredFitError < 20), 0:0.1:10);
        bar(xout,n)
        hold on
        [n,xout] = hist(fitError(fitError > 0 & fitError < 20), 0:0.1:10);
        bar(xout,n, 'green')
%         [n,xout] = hist(badFitError(badFitError > 0 & badFitError < 20), 0:0.1:10);
%         bar(xout,n, 'red')
        hold off
        xlabel('Fit Error'); ylabel('Frequency');
        title('Fit Error');
        legend('unfiltered','initially included')
        xlim([0 8]);
        
        subplot(2,2,3)
        hist(catPSFfits(initGoodFits,ampRatioCol), 100)
        [n,xout] = hist(catPSFfits(:,ampRatioCol), 0:0.01:1);
        bar(xout,n)
        hold on
        [n,xout] = hist(catPSFfits(initGoodFits,ampRatioCol), 0:0.01:1);
        bar(xout,n, 'green')
        hold off
        xlabel('Amplitude Ratio'); ylabel('Frequency');
        title('Amplitude Ratio');
        xlim([-0.1 1]);
        legend('unfiltered','initially included')
        
        subplot(2,2,4)
        hist(catPSFfits(initGoodFits,sigmaRatioCol), 100)
        xlabel('Sigma Ratio'); ylabel('Frequency');
        title('Sigma Ratio');
        xlim([-0.1 1]);

    end
    %% Chose a desired parameter set for reconstruction
    
    dlg_title = 'Please Input Parameters';
    prompt = {  'Size of points in reconstruction',...
        'Temporal Color Coding',...
        'White light shift X (in nm)',...
        'White light shift Y (in nm)',...
        'Z range lower bound(in nm)',...
        'Z range upper bound(in nm)',...
        'First frame',...
        'Last frame'...
        'Plot as Tracks'...
        'Display templates in points'...
        };
    def = { ...
        num2str(scatterSize), ...
        num2str(useTimeColors), ...
        num2str(wlShiftX), ...
        num2str(wlShiftY), ...
        num2str(zRange(1)), ...
        num2str(zRange(2)), ...
        num2str(frameRange(1)), ...
        num2str(frameRange(2)), ...
        num2str(plotAsTracks), ...
        num2str(showTemplates)...
        };
    num_lines = 1;
    inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
        
    scatterSize = str2double(inputdialog{1});
    useTimeColors = str2double(inputdialog{2});
    wlShiftX = str2double(inputdialog{3});
    wlShiftY = str2double(inputdialog{4});
    zRange = [str2double(inputdialog{5}) str2double(inputdialog{6})];
    frameRange = [str2double(inputdialog{7}) str2double(inputdialog{8})];
    plotAsTracks = str2double(inputdialog{9});
    showTempates = str2double(inputdialog{10});
    dlg_title = 'Please Input Parameters';
    prompt = {  ...
%         'Lobe sigma lower bound (in pixel)',...
%         'Lobe sigma upper bound (in pixel)',...
        'Lobe distance lower bound (in pixel)',...
        'Lobe distance upper bound (in pixel)',...
        'Amplitude ratio limit',...
        'Sigma ratio limit',...
        'Photon weighted fit error lower bound',...
        'Photon weighted fit error upper bound',...
        'Number of photons lower bound',...
        'Number of photons upper bound',...
        'xyPrec lower bound (nm)',...
        'xyPrec upper bound (nm)',...
        'zPrec lower bound (nm)',...
        'zPrec upper bound (nm)',...
        'Minimum template match values',...
        };
    def = { ...
%         num2str(sigmaBounds(1)), ...
%         num2str(sigmaBounds(2)), ...
        num2str(lobeDistBounds(1)), ...
        num2str(lobeDistBounds(2)), ...
        num2str(ampRatioLimit), ...
        num2str(sigmaRatioLimit), ...
        num2str(fitErrorRange(1)), ...
        num2str(fitErrorRange(2)), ...
        num2str(numPhotonRange(1)), ...
        num2str(numPhotonRange(2)), ...
        num2str(xyPrecRange(1)),...
        num2str(xyPrecRange(2)),...
        num2str(zPrecRange(1)),...
        num2str(zPrecRange(2)),...
        ['[' num2str(threshVals) ']'],...
        };
    num_lines = 1;
    inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
    
    %     sigmaBounds = [str2double(inputdialog{1}) str2double(inputdialog{2})];   %[1.2 2.7];    % sets [min max] allowed sigma for double Gaussian fit (in units of pixels)
    lobeDistBounds = [str2double(inputdialog{1}) str2double(inputdialog{2})];  %[7.0 9.5]; % sets [min max] allowed interlobe distance for double Gaussian fit (in units of pixels)
    ampRatioLimit = str2double(inputdialog{3});
    sigmaRatioLimit = str2double(inputdialog{4});
    fitErrorRange = [str2double(inputdialog{5}) str2double(inputdialog{6})];
    numPhotonRange = [str2double(inputdialog{7}) str2double(inputdialog{8})];
    xyPrecRange = [str2double(inputdialog{9}) str2double(inputdialog{10})];
    zPrecRange = [str2double(inputdialog{11}) str2double(inputdialog{12})];
    threshVals = str2num(inputdialog{13});
    %% Plot the white light image if specified
    if any(whiteLightFile ~= 0)
        if ~exist('whiteLight')
        whiteLightInfo = imfinfo(whiteLightFile);
        whiteLight = zeros(whiteLightInfo(1).Height, whiteLightInfo(1).Width);
        % average white light images together to get better SNR
        for a = 1:length(whiteLightInfo)
            whiteLight = whiteLight + double(imread(whiteLightFile,a, ...
                'Info', whiteLightInfo));
        end
        % resize white light to the size of the ROI of the single molecule fits
        whiteLight = whiteLight(ROI_initial(2):ROI_initial(2)+ROI_initial(4)-1,ROI_initial(1):ROI_initial(1)+ROI_initial(3)-1);
        % rescale white light image to vary from 0 to 1
        %         whiteLight = (whiteLight-min(whiteLight(:)))/(max(whiteLight(:))-min(whiteLight(:)));
        borderStart = 80; borderEnd = 230;
            wlCenter = whiteLight(borderStart:borderEnd,...
                borderStart:borderEnd);
            whiteLight = (whiteLight-min(wlCenter(:)))/...
                (max(wlCenter(:))-min(wlCenter(:)));
            whiteLight(whiteLight(:)>1) = 1; whiteLight(whiteLight(:)<0) = 0;
            
        end
        [xWL, yWL] = meshgrid((ROI_initial(1):ROI_initial(1)+ROI_initial(3)-1) * nmPerPixel + wlShiftX, ...
            (ROI_initial(2):ROI_initial(2)+ROI_initial(4)-1) * nmPerPixel + wlShiftY);
        %         [xWL yWL] = meshgrid((1:(whiteLightInfo(1).Width)) * nmPerPixel + wlShiftX, ...
        %         (1:(whiteLightInfo(1).Height)) * nmPerPixel + wlShiftY);
    end
    

    
    %% Now re-evaluate the goodness of the fits
    
    % Conditions for fits (play with these):
    % (1) Amplitude of both lobes > 0
    % (2) All locations x1,y1, x2,y2 lie inside area of small box
    % (3) All sigmas need to be > sigmaBound(1) and < sigmaBound(2)
    % (4) Distance between lobes needs to be > lobeDist(1) pixels and < lobeDist(2) pixels
    % (5) Make sure amplitudes are within 100% of one another
    % (6) Make sure totalFitError/(total number of photons) is within the fitErrorRange
    
    fitErrorCol = 16;
    goodFitFlagCol = 17;
    numPhotonCol = 21;
    bkgndCol = 15;
    lobeDistCol = 22;
    ampRatioCol = 23;
    sigmaRatioCol = 24;
        
    for i = 1:size(catPSFfits,1)
        
        % compute localization precision as a function of the number of photons
        % Empirically determined amplitudes for fitting function based on
        % localization precisision calibration collected on 20120518 on 8a back setup.
        amplitude =  [  361035.867260138,22.2956414971275;...   %   [A1x  A2x]
            348907.934759022,28.3183226442783;...   %   [A1y  A2y]
            840446.405407229,23.3314294806927];      %   [A1z  A2z]
        
        numPhotons = catPSFfits(i,numPhotonCol);
        meanBkgnd = catPSFfits(i,bkgndCol);
        
        if any(meanBkgnd<0)
            meanBkgnd(meanBkgnd<0) = 0;
        end
        % Equation 4 of Stallinga and Rieger, ISBI, Barcelona conference proveedings
        sigmaX = sqrt(amplitude(1,1) .* (1./numPhotons) + amplitude(1,1)*4*amplitude(1,2) .* meanBkgnd./(numPhotons).^2 + amplitude(1,1) .* (1./numPhotons) .* sqrt((2*amplitude(1,2)*(meanBkgnd./numPhotons))./(1+(4*amplitude(1,2)*(meanBkgnd./numPhotons)))));
        sigmaY = sqrt(amplitude(2,1) .* (1./numPhotons) + amplitude(2,1)*4*amplitude(2,2) .* meanBkgnd./(numPhotons).^2 + amplitude(2,1) .* (1./numPhotons) .* sqrt((2*amplitude(2,2)*(meanBkgnd./numPhotons))./(1+(4*amplitude(2,2)*(meanBkgnd./numPhotons)))));
        sigmaZ = sqrt(amplitude(3,1) .* (1./numPhotons) + amplitude(3,1)*4*amplitude(3,2) .* meanBkgnd./(numPhotons).^2 + amplitude(3,1) .* (1./numPhotons) .* sqrt((2*amplitude(3,2)*(meanBkgnd./numPhotons))./(1+(4*amplitude(3,2)*(meanBkgnd./numPhotons)))));
        
        
        
        if catPSFfits(i,goodFitFlagCol) == -1001 || ...
                catPSFfits(i,goodFitFlagCol) == -1002 || ...
                catPSFfits(i,goodFitFlagCol) == -1003
            
            continue
            
        %gaussian sigma ratio filter
        elseif catPSFfits(i,sigmaRatioCol) > sigmaRatioLimit;
            
            catPSFfits(i,goodFitFlagCol) = -1004;
        
        % lobe distance filter
        elseif catPSFfits(i,lobeDistCol) < lobeDistBounds(1) || catPSFfits(i,lobeDistCol) > lobeDistBounds(2)
            
            catPSFfits(i,goodFitFlagCol) = -1005;
       
        % amplitude ratio filter
        elseif catPSFfits(i,ampRatioCol) > ampRatioLimit;
            
            catPSFfits(i,goodFitFlagCol) = -1006;
        
        % weighted error filter
        elseif catPSFfits(i,fitErrorCol)*conversionFactor/catPSFfits(i,numPhotonCol) > fitErrorRange(2)  || ...
                catPSFfits(i,fitErrorCol)*conversionFactor/catPSFfits(i,numPhotonCol) < fitErrorRange(1)
            
            catPSFfits(i,goodFitFlagCol) = -1007;
        
        % localization precision filter     
        elseif sigmaX < xyPrecRange(1) || sigmaX > xyPrecRange(2) ||...
               sigmaZ < zPrecRange(1) || sigmaZ > zPrecRange(2)
            
            catPSFfits(i,goodFitFlagCol) = -1008;
        
        % template match strength filter
        elseif catPSFfits(i,templateStrCol) < threshVals(catPSFfits(i,templateNumCol))
            catPSFfits(i,goodFitFlagCol) = -1009;
            
        else
            catPSFfits(i,goodFitFlagCol) = 3;
        end
        
    end
    
    %% load valid xyz locations
    
%     goodFits = false(size(totalPSFfits,1),1);
    
    goodFits = catPSFfits(:,17) > 0; % totalPSFfits(:,17) > -inf;
    goodFits = goodFits & ~isnan(catPSFfits(:,25));%%%%
    badFits = catPSFfits(:,17) < 0;
    goodFits = goodFits & catPSFfits(:,1) >= frameRange(1) & catPSFfits(:,1) <= frameRange(2);
    goodFitsNoPhotFilt = sum(goodFits);
    goodFits = goodFits & (catPSFfits(:,numPhotonCol) >= numPhotonRange(1)) & (catPSFfits(:,numPhotonCol) <= numPhotonRange(2));
    if goodFitsNoPhotFilt - sum(goodFits) > goodFitsNoPhotFilt / 2
        warning(['More than half of the fits were thrown out due to '...
                'restrictions on the # photons. Double check this limit.']);
    end
    if sum(goodFits) < 5
        warning('Very few (<5) fits passed the filters. Double-check limits.');
    end
    % corrects zRange for index mismatch (see below for the inverse
    % transformation to the z position)
    corrzRange = zRange * nOil/nSample;
    
    % extract data; check if fiducial correction is in play (i.e., col 30 is not all nans)
    % notate fid-corrected localizations separately
    if sum(isnan(catPSFfits(:,30))) ~= size(catPSFfits,1)
        goodFits = goodFits & catPSFfits(:,30) >= corrzRange(1) & catPSFfits(:,30) <= corrzRange(2); % try to fid-corrected values for z range
        xLocPix = catPSFfits(goodFits,18)/nmPerPixel;
        yLocPix = catPSFfits(goodFits,19)/nmPerPixel;
        xLoc = catPSFfits(goodFits,28);
        yLoc = catPSFfits(goodFits,29);
        zLoc = catPSFfits(goodFits,30);
        tempNumLoc = catPSFfits(goodFits,templateNumCol);
        
        xLoc_bad = catPSFfits(badFits,28);
        yLoc_bad = catPSFfits(badFits,29);
        % still generate xLoc etc. as below, but call them 'raw' if fids
        % were used. these can then be used for registration, where both
        % focal shift correction and fiducial correction should be
        % performed AFTER registering the two channels.
        xLocPixRaw = catPSFfits(goodFits,18)/nmPerPixel;
        yLocPixRaw = catPSFfits(goodFits,19)/nmPerPixel;
        xLocRaw = catPSFfits(goodFits,25);
        yLocRaw = catPSFfits(goodFits,26);
        zLocRaw = catPSFfits(goodFits,27);
        xLoc_badRaw = catPSFfits(badFits,25);
        yLoc_badRaw = catPSFfits(badFits,26);
        
    else % if no fid-corrected traces were generated (column 30 is nan)
        goodFits = goodFits & catPSFfits(:,27) >= corrzRange(1) & catPSFfits(:,27) <= corrzRange(2);        
        xLocPix = catPSFfits(goodFits,18)/nmPerPixel;
        yLocPix = catPSFfits(goodFits,19)/nmPerPixel;
        xLoc = catPSFfits(goodFits,25);
        yLoc = catPSFfits(goodFits,26);
        zLoc = catPSFfits(goodFits,27);
        xLoc_bad = catPSFfits(badFits,25);
        yLoc_bad = catPSFfits(badFits,26);
        tempNumLoc = catPSFfits(goodFits,templateNumCol);
        
    end
        
    clear corrzRange
    zLoc_IndexCorrected = zLoc * nSample/nOil;
    if exist('xLocRaw','var')
        zLoc_IndexCorrectedRaw = zLocRaw * nSample/nOil;
    end
    
    numPhotons = catPSFfits(goodFits,21);
    lobeDist = catPSFfits(goodFits,22);
    ampRatio = catPSFfits(goodFits,23);
    sigmaRatio = catPSFfits(goodFits,24);
    
%     meanBkgnd = totalPSFfits(goodFits,15)*conversionFactor;
    meanBkgnd = catPSFfits(goodFits,15);  % output from template match is already in units of photons.
    frameNum = catPSFfits(goodFits,1);
    PSFfits_bad = catPSFfits(badFits,:);
    
    %% ask user what region to plot in superresolution image
    
    h2Dfig=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
    if whiteLightFile ~= 0
        xRange = xWL(1,:);
        yRange = yWL(:,1);
        % pick region that contains background
        imagesc(xRange,yRange,whiteLight, [0 1]);axis image;colormap gray;hold on;
    end
    
    % plot is faster than scatter
    if showTemplates
        grayBitSize = 26;
        templateColors = [gray(grayBitSize); lines(6)];
        for blah = 1:6
            toPlot = tempNumLoc == blah;
        plot(xLoc(toPlot),yLoc(toPlot),'.','MarkerSize',5,'Color',templateColors(blah+grayBitSize,:))
        end
        cbar_handle = colorbar('location','eastoutside');
        colormap(templateColors);
        caxis([0.5 1.1])
        set(cbar_handle,'YTickLabel',{'1','2','3','4','5','6'},...%num2str(max(Z0(:)))},...
        'ytick',linspace(1.1-0.6*6/32+0.6*6/32/12,1.1-0.6*6/32/12,6),'ylim',[1.1-0.6*6/32 1.1]);
        title(cbar_handle,'Template Number')
    else
    plot(xLoc,yLoc,'.','MarkerSize',1);
    end
    xlim([min(xLoc)-500, max(xLoc)+500]);
    ylim([min(yLoc)-500, max(yLoc)+500]);
    xlabel('x (nm)');ylabel('y (nm)');
    axis ij;
    axis square;
    
    if pass == 1
        ROI = imrect(gca,[min(xLoc(:)) min(yLoc(:)) max(xLoc(:))-min(xLoc(:)) max(yLoc(:))-min(yLoc(:))]);
    else
        ROI = imrect(gca,[ROI(1) ROI(2) ROI(3) ROI(4)]);
    end
    
    title({'Double-click to choose region that will be plotted in 3D scatterplot' ...
        mat2str(ROI.getPosition)});
    addNewPositionCallback(ROI,@(p) title({'Double-click to choose region that will be plotted in 3D scatterplot' ...
        ['[xmin ymin width height] = ' mat2str(p,3)]}));
    % make sure rectangle stays within image bounds
    fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
    setPositionConstraintFcn(ROI,fcn);
    ROI = wait(ROI);
    clear avgImg fcn
    
    %% filter out localizations outside of ROI
    
    validPoints = xLoc>=ROI(1) & xLoc<=ROI(1)+ROI(3) & yLoc>ROI(2) & yLoc<=ROI(2)+ROI(4) & numPhotons>0;
    invalidPoints = xLoc_bad>=ROI(1) & xLoc_bad<=ROI(1)+ROI(3) & yLoc_bad>ROI(2) & yLoc_bad<=ROI(2)+ROI(4) ;
    
    if ~any(validPoints)
        disp('You chose an area without any points');
        continue
    end
    
    xLocPix = xLocPix(validPoints);
    yLocPix = yLocPix(validPoints); 
    xLoc = xLoc(validPoints);
    yLoc = yLoc(validPoints);
    zLoc = zLoc(validPoints);
    zLoc_IndexCorrected = zLoc_IndexCorrected(validPoints);
    
    if exist('xLocRaw','var')
        xLocPixRaw = xLocPixRaw(validPoints);
        yLocPixRaw = yLocPixRaw(validPoints); 
        xLocRaw = xLocRaw(validPoints);
        yLocRaw = yLocRaw(validPoints);
        zLocRaw = zLocRaw(validPoints);
        zLoc_IndexCorrectedRaw = zLoc_IndexCorrectedRaw(validPoints);
    end
    
%     [std(xLoc) std(yLoc) std(zLoc)]
    numPhotons = numPhotons(validPoints);
    meanBkgnd = meanBkgnd(validPoints);
    frameNum = frameNum(validPoints);
    lobeDist = lobeDist(validPoints);
    ampRatio = ampRatio(validPoints);
    sigmaRatio = sigmaRatio(validPoints);
    
    PSFfits_bad = PSFfits_bad(invalidPoints,:);
    hRejections = figure;
    subplot(2,2,1:2)
    x = -1008:1:-1001;
    hist(PSFfits_bad(PSFfits_bad(:,17)<-10,17),x)
    xlabel('Error Flag');ylabel('Frequency');
    title({[num2str(size(PSFfits_bad,1)) ' bad localizations']});
    
    %% compute localization precision as a function of the number of photons
    % ToDo:  Repeat this calibration
    
    % Empirically determined amplitudes for fitting function based on
    % localization precisision calibration in Nano Letters Paper.
    %     sigmaX = 410./(1.5*numPhotons./sqrt(meanBkgnd)).^0.47;
    %     sigmaY = 550./(1.5*numPhotons./sqrt(meanBkgnd)).^0.52;
    %     sigmaZ = 829./(1.5*numPhotons./sqrt(meanBkgnd)).^0.49;
    
    %     % Empirically determined amplitudes for fitting function based on
    %     % localization precisision calibration collected on 20120402 on 8a back setup.
    %     amplitude =  [  606316.910875840,1351845.90313904;...   %   [A1x  A2x]
    %                     463419.307260597,1230505.00679917;...   %   [A1y  A2y]
    %                     990499.159483260,3178237.19926875]      %   [A1z  A2z]
    
    % Empirically determined amplitudes for fitting function based on
    %     % localization precisision calibration collected on 20120518 on 8a back setup.
    %     %%
    %     numPhotons = 500;
    %     meanBkgnd = 5;
    %     amplitude =  [  467376.158647402,1696621.37143132;...   %   [A1x  A2x]
    %                     472618.572573088,1601434.40940051;...   %   [A1y  A2y]
    %                     1096609.27949884,3959185.51690004];      %   [A1z  A2z]
    %
    %     sigmaX = sqrt(amplitude(1,1) .* (1./numPhotons).^1 + amplitude(1,2) .* (meanBkgnd./numPhotons).^2)
    %     sigmaY = sqrt(amplitude(2,1) .* (1./numPhotons).^1 + amplitude(2,2) .* (meanBkgnd./numPhotons).^2)
    %     sigmaZ = sqrt(amplitude(3,1) .* (1./numPhotons).^1 + amplitude(3,2) .* (meanBkgnd./numPhotons).^2)
    %     %%
    
    amplitude =  [  361035.867260138,22.2956414971275;...   %   [A1x  A2x]
        348907.934759022,28.3183226442783;...   %   [A1y  A2y]
        840446.405407229,23.3314294806927];      %   [A1z  A2z]
    
    % Equation 4 of Stallinga and Rieger, ISBI, Barcelona conference proveedings
    if any(meanBkgnd<0);
        correctedBG = find(meanBkgnd<0);
        warning([num2str(length(correctedBG))...
            'values of meanBkgnd (indices below) were negative! Changing to 0.']);
        meanBkgnd(correctedBG)=0;
    end
    
    sigmaX = sqrt(amplitude(1,1) .* (1./numPhotons) + amplitude(1,1)*4*amplitude(1,2) .* meanBkgnd./(numPhotons).^2 + amplitude(1,1) .* (1./numPhotons) .* sqrt((2*amplitude(1,2)*(meanBkgnd./numPhotons))./(1+(4*amplitude(1,2)*(meanBkgnd./numPhotons)))));
    sigmaY = sqrt(amplitude(2,1) .* (1./numPhotons) + amplitude(2,1)*4*amplitude(2,2) .* meanBkgnd./(numPhotons).^2 + amplitude(2,1) .* (1./numPhotons) .* sqrt((2*amplitude(2,2)*(meanBkgnd./numPhotons))./(1+(4*amplitude(2,2)*(meanBkgnd./numPhotons)))));
    sigmaZ = sqrt(amplitude(3,1) .* (1./numPhotons) + amplitude(3,1)*4*amplitude(3,2) .* meanBkgnd./(numPhotons).^2 + amplitude(3,1) .* (1./numPhotons) .* sqrt((2*amplitude(3,2)*(meanBkgnd./numPhotons))./(1+(4*amplitude(3,2)*(meanBkgnd./numPhotons)))));
    
    
    %     %%
    %     amplitude =  [  360000,22;...   %   [A1x  A2x]
    %                     350000,28;...   %   [A1y  A2y]
    %                     840000,23];      %   [A1z  A2z]
    %     numPhotons = 4000;
    %     meanBkgnd = 7;
    %     sigmaX - sqrt(amplitude(1,1) .* (1./numPhotons) + amplitude(1,1)*4*amplitude(1,2) .* meanBkgnd./(numPhotons).^2 + amplitude(1,1) .* (1./numPhotons) .* sqrt((2*amplitude(1,2)*(meanBkgnd./numPhotons))./(1+(4*amplitude(1,2)*(meanBkgnd./numPhotons)))))
    %     sigmaY - sqrt(amplitude(2,1) .* (1./numPhotons) + amplitude(2,1)*4*amplitude(2,2) .* meanBkgnd./(numPhotons).^2 + amplitude(2,1) .* (1./numPhotons) .* sqrt((2*amplitude(2,2)*(meanBkgnd./numPhotons))./(1+(4*amplitude(2,2)*(meanBkgnd./numPhotons)))))
    %     sigmaZ - sqrt(amplitude(3,1) .* (1./numPhotons) + amplitude(3,1)*4*amplitude(3,2) .* meanBkgnd./(numPhotons).^2 + amplitude(3,1) .* (1./numPhotons) .* sqrt((2*amplitude(3,2)*(meanBkgnd./numPhotons))./(1+(4*amplitude(3,2)*(meanBkgnd./numPhotons)))))
    %     %%
    
    
    meanNumPhotons = mean(numPhotons);
    %     localizationPrecision = [mean(sigmaX),mean(sigmaY),mean(sigmaZ)];
    %     frameNum;
    subplot(2,2,3:4)
    %     hist(frameNum,1:length(frames))
    hist(frameNum,frameNum(1):frameNum(length(frameNum)))
    ylim([0 1.4])
    xlabel('Frame Number');ylabel('Single Molecule Fit')
    
    %% plot statistics on these localizations
    
    hStatsFig=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
    subplot(1,4,1);
    hist(numPhotons,round(length(xLoc)/20));
    xlabel('Number of photons per localization');
    subplot(1,4,2);
    hist(sigmaX,round(length(xLoc)/20));
    xlabel('\sigma_x (nm)');
    subplot(1,4,3);
    hist(sigmaY,round(length(xLoc)/20));
    xlabel('\sigma_y (nm)');
    subplot(1,4,4);
    hist(sigmaZ,round(length(xLoc)/20));
    xlabel('\sigma_z (nm)');
    
    %imwrite(frame2im(getframe(h)),[filePrefix ' localization stats.tif']);
    
    %% plot 3D scatterplot of localizations with white light
    
    h3Dfig = figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','k','renderer','opengl', 'Toolbar', 'figure'); hold on;
    if whiteLightFile~=0
        %imagesc(xRange,yRange,whiteLight);axis image;colormap gray;hold on;
        [x,y,z] = meshgrid(xRange,yRange,[min(zLoc_IndexCorrected) max(zLoc_IndexCorrected)]);
        xslice = []; yslice = []; zslice = min(zLoc_IndexCorrected);
        h=slice(x,y,z,repmat(whiteLight,[1 1 2]),xslice,yslice,zslice,'nearest');
        set(h,'EdgeColor','none','FaceAlpha',0.75);
        colormap gray;
    end
    
    if plotAsTracks == 1
        markerColors = jet(frameNum(length(frameNum))-frameNum(1)+1);
        
        for a = 1:length(frameNum)-1
        
%             plot3(xLoc,yLoc,zLoc_IndexCorrected,'-',...
%                 'Color',[1 1 0]);
            plot3(xLoc(a:a+1),yLoc(a:a+1),zLoc_IndexCorrected(a:a+1),'LineWidth',scatterSize/12,...
                'Color',markerColors(frameNum(a)-frameNum(1)+1,:));
            
%             scatter3(xLoc(a),yLoc(a),zLoc_IndexCorrected(a),scatterSize,'filled',...
%                 'MarkerFaceColor', markerColors(frameNum(a)-frameNum(1)+1,:),...
%                 'MarkerEdgeColor', markerColors(frameNum(a)-frameNum(1)+1,:));
        end
    elseif useTimeColors == 0
        % plot is faster than scatter
        plot3(xLoc,yLoc,zLoc_IndexCorrected,'.','MarkerSize',scatterSize/3,...
            'Color',[1 1 0]);
    %     scatter3(xLoc,yLoc,zLoc,scatterSize,[1 1 0],'filled');
    elseif useTimeColors == 1
        %         scatter3(xLoc(a),yLoc(a),zLoc(a),scatterSize,frameNum(1):frameNum(length(frameNum)),'filled')
        markerColors = jet(frameNum(length(frameNum))-frameNum(1)+1);
        %     for a = 1:length(frameNum)
        %         scatter3(xLoc(a)-min(xLoc(:)),yLoc(a)-min(yLoc(:)),zLoc(a),scatterSize,'filled',...
        %             'MarkerFaceColor', markerColors(frameNum(a)-frameNum(1)+1,:),...
        %             'MarkerEdgeColor', markerColors(frameNum(a)-frameNum(1)+1,:));
        %     end
        
        for a = 1:length(frameNum)
            scatter3(xLoc(a),yLoc(a),zLoc_IndexCorrected(a),scatterSize,'filled',...
                'MarkerFaceColor', markerColors(frameNum(a)-frameNum(1)+1,:),...
                'MarkerEdgeColor', markerColors(frameNum(a)-frameNum(1)+1,:));
        end
    end
    axis vis3d equal;
    %     xlim([min(xLoc(:)) max(xLoc(:))]);
    %     ylim([min(yLoc(:)) max(yLoc(:))]);
    %     xlim([min(xLoc(:)) max(xLoc(:))]-min(xLoc(:)));
    %     ylim([min(yLoc(:)) max(yLoc(:))]-min(yLoc(:)));
    xlim([min(xLoc(:)) max(xLoc(:))]);
    ylim([min(yLoc(:)) max(yLoc(:))]);
    xlabel('x (nm)');ylabel('y (nm)');zlabel('z (nm)');
    
    ROICenterX = ROI(1)+ROI(3)/2;
    ROICenterY = ROI(2)+ROI(4)/2;
    
    % todo: add a way to get around this if
    if exist('laser_x_nm', 'var')
        
        distToPeakIntensity = sqrt((laser_x_nm-ROICenterX)^2 + (laser_y_nm-ROICenterY)^2);
        meanIntensityInROI = peakIntensity * exp(-((2*distToPeakIntensity^2)/((2*mean([sigma_x_nm, sigma_y_nm]))^2)));
        
        title({[num2str(length(xLoc)) ' localizations'];...
            ['Mean Number of Signal Photons = ' num2str(meanNumPhotons) ' per frame'];...
            ['Mean Number of Background Photons = ' num2str(mean(meanBkgnd)) ' per pixel per frame'];...
            ['Localization Precision \sigma_x = ' num2str(mean(sigmaX)) ' nm'];...
            ['Localization Precision \sigma_y = ' num2str(mean(sigmaY)) ' nm'];...
            ['Localization Precision \sigma_z = ' num2str(mean(sigmaZ)) ' nm'];...
            ['Laser Intensity = ' num2str(meanIntensityInROI) ' W/cm^2']},...
            'color','w');
    else
        
        title({[num2str(length(xLoc)) ' localizations'];...
            ['Mean Number of Signal Photons = ' num2str(meanNumPhotons) ' per frame'];...
            ['Mean Number of Background Photons = ' num2str(mean(meanBkgnd)) ' per pixel per frame'];...
            ['Localization Precision \sigma_x = ' num2str(mean(sigmaX)) ' nm'];...
            ['Localization Precision \sigma_y = ' num2str(mean(sigmaY)) ' nm'];...
            ['Localization Precision \sigma_z = ' num2str(mean(sigmaZ)) ' nm']},...
            'color','w');
    end

    set(gca,'color','k');
    set(gca,'xcolor','w');set(gca,'ycolor','w');set(gca,'zcolor','w');
    
    %% Construct a questdlg with three options
    
    % f = figure;
    h = uicontrol('Position',[20 20 200 40],'String','Continue',...
        'Callback','uiresume(gcbf)');
    % disp('This will print immediately');
    uiwait(gcf);
    % disp('This will print after you click Continue');
%     close(f);
    
    dlg_title = 'Replot';
    prompt = {'Would you like to replot with a different parameter set?'};
    def =       { 'Yes'  };
    questiondialog = questdlg(prompt,dlg_title, def);
    % Handle response
    switch questiondialog
        case 'Yes'
            pass = pass + 1;
            close 
            close
            close
            close 
        case 'No'
            anotherpass = false;
        case 'Cancel'
            error('User cancelled the program');
    end
        
end

%% prompt to save data
[saveFile, savePath] = uiputfile({'*.*'},'Enter a directory title for this ROI. Otherwise, click cancel.');
savePath = [savePath saveFile filesep];
mkdir(savePath);
if ~isequal(saveFile,0)
    save([savePath 'Output'],'xLocPix','yLocPix','xLoc','yLoc','zLoc','zLoc_IndexCorrected','numPhotons','meanBkgnd','sigmaX','sigmaY','sigmaZ','frameNum',...
        'zRange','frameRange','sigmaBounds','lobeDistBounds','ampRatioLimit','sigmaRatioLimit','fitErrorRange','numPhotonRange',...
        'lobeDist','ampRatio','sigmaRatio','wlShiftX', 'wlShiftY','goodFits','fidTrackX', 'fidTrackY', 'fidTrackZ', 'nmPerPixel','whiteLightFile','threshVals',...
        'spatialCorr','useCurrent','currFidIdx','catPSFfits','fidsToUse','fidChannel','imgShiftX','imgShiftY','ROI');
    if exist('xLocRaw');
        save([savePath 'Output'],'xLocRaw','yLocRaw','zLocRaw','zLoc_IndexCorrectedRaw','-append');
    end
    if exist('whiteLight');
        save([savePath 'Output'],'whiteLight','xWL','yWL','-append');
    end
    if exist('sigma_x_nm');
        save([savePath 'Output'],'powerAtObjective','sigma_x_nm','sigma_y_nm','laser_x_nm','laser_y_nm','theta','peakIntensity','-append');
    end
end
%%
% output excel spreadsheet
% textHeader = {'frame number' ...
%     'raw x location (pix)' ...
%     'raw y location (pix)' ...    
%     'fiduciary corrected x location (nm)' ...
%     'fiduciary corrected y location (nm)' ...
%     'fiduciary corrected z location (nm)' ...
%     'sigma x (nm)' ...
%     'sigma y (nm)' ...
%     'sigma z (nm)' ...
%     'number of photons' ...
%     'mean background photons' };
% output = [frameNum, xLocPix yLocPix xLoc, yLoc, zLoc, sigmaX, sigmaY, sigmaZ, numPhotons, meanBkgnd];
% xlswrite([savePath saveFile(1:length(saveFile)-4) '.xlsx'], [textHeader; ...
%     num2cell(output)], ...
%     'valid PSF fits');

% Save figures  

saveas(h3Dfig,[savePath '3D.fig']);
close(h3Dfig)
saveas(hStatsFig,[savePath 'stats.fig']);
close(hStatsFig)
saveas(hRejections,[savePath 'rejections.fig']);
close(hRejections)
saveas(h2Dfig,[savePath '2D.fig']);
close(h2Dfig)

% print(gcf,'-depsc','-r2400','-loose',[savePath saveFile(1:length(saveFile)-4) '_2D.eps']);

end
