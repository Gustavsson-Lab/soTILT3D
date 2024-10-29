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

function [totalPSFfits, numFrames, fidTrackX, fidTrackY, fidTrackZ,spatialCorr,useCurrent,fidsToUse,fidChannel,imgShiftX,imgShiftY] = ...
    f_concatSMfits(fitFilePrefix,useFidCorrections,fidFilePrefix,logFile,logPath,channel,calFile,currFidIdx,nmPerPixel,smacmFile)
%clear all;
% close all;
numSyncFrames = 25;
useDenoising = 1;
undoZangleZero = 1;
%Construct a questdlg with three options
                dlg_title = 'Use local cals?';
                prompt = { 'Do you want to try to use spatially dependent calibrations?' };
                def =       { 'No'  };
                questiondialog = questdlg(prompt,dlg_title, def);
                % Handle response
                switch questiondialog
                    case 'Yes'
                        spatialCorr=true;
                    case 'No'
                        spatialCorr=false;
                    case 'Cancel'
                        error('User cancelled the program');
                end

                dlg_title = 'Correct for image offset?';
    prompt = {...
        'X offset of image (starting pixel)',...
        'Y offset of image (starting pixel)'...
             };
    def = {'1','1'}; % e.g. if 320x320 of 512x512 in far corner, use 193, 193
    num_lines = 1;
    inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
    imgShiftX = str2double(inputdialog{1})-1;
    imgShiftY = str2double(inputdialog{2})-1;
    
    
    
% dlg_title = 'Flatten field?';
% prompt = { 'Do you want to subtract out angle to make z = 0 when theta = 0?' };
% def =       { 'Yes'  };
% questiondialog = questdlg(prompt,dlg_title, def);
% % Handle response
% switch questiondialog
%     case 'Yes'
%         undoZangleZero=false;
%     case 'No'
%         undoZangleZero=true;
%     case 'Cancel'
%         error('User cancelled the program');
% end
undoZangleZero = true;
disp('Calibrations are being reset, so z = 0 is not necessarily theta = 0.')
% Can use sif logs to filter which frames are included for fiducial correction.
% If logPath set to 0, all frames are used (every frame with a visible fiducial,
% even when the 'wrong' laser is on).
% Set to 0 because it is not clear that using filtering is ideal: sometimes 
% using the frames with the 'wrong laser' can introduce a slight bias, but 
% sometimes (especially for z traces) having extra frames is worth it,
% unless we can find a good way for waveletFidTracks to 'estimate' what's
% going on in large gaps

logFile=0; logPath=0; %channel=[];

%% prepare fiducial data / blank out fiducial arrays
if useFidCorrections
    %% load raw fiduciary data
    % [fidFile fidPath] = uigetfile({'*.mat';'*.*'},'Open data file #1 with raw fiduciary fits');
    for fileNum=1:length(fidFilePrefix)
        %         fidFiles = [fidFiles; {[fidPath fidFile]}];
        % load data
        load([fidFilePrefix{fileNum} 'raw fits.mat'],'PSFfits','numFrames','numMoles');
        
        if fileNum == 1
            tempPSFfits = PSFfits(:,1:23);
            numFramesInFiles = numFrames;
        else
            numFramesInFiles = [numFramesInFiles numFrames];
            PSFfits(:,1) = PSFfits(:,1) + sum(numFramesInFiles(1:fileNum-1));
            tempPSFfits = [tempPSFfits; PSFfits(:,1:23)];
        end
        
        %         fileNum = fileNum+1;
        %         [fidFile fidPath] = uigetfile({'*.mat';'*.*'},...
        %             ['Open data file #' num2str(fileNum) ' with raw fiduciary fits']);
    end
    
    PSFfits = shiftData(tempPSFfits,'FID',imgShiftX,imgShiftY,nmPerPixel);
    numFrames = sum(numFramesInFiles);
%     clear tempPSFfits;
        
    %% prep data file for use in picking fiducials
    dataFileInfo = imfinfo(smacmFile);
    dataWidth = dataFileInfo(1).Width;
    dataHeight = dataFileInfo(1).Height;
    dataMockup = zeros(dataHeight,dataWidth);
    
    dataLength = length(dataFileInfo);
    goodFrames = unique(PSFfits(PSFfits(:,13)>0,1));
    
    for a = goodFrames(1:10)' % output averaged data to pick fids
        dataMockup = dataMockup + double(imread(smacmFile,a,'Info',dataFileInfo));
    end
    dataMockup = dataMockup/10;
    
    %% select fiducial(s) to use
    
    xCol = 14; yCol = 15;
    moleLocs = nan(numMoles,2);
    for j = 1:numMoles
        moleLocs(j,1) = nanmean(PSFfits(PSFfits(:,2)==j,xCol));
        moleLocs(j,2) = nanmean(PSFfits(PSFfits(:,2)==j,yCol));
    end
    moleLocs = moleLocs / nmPerPixel;
    
    xRange = imgShiftX + [1 dataWidth];
    yRange = imgShiftY + [1 dataHeight];
    hSelFig = figure; hold on;
    imagesc(yRange,xRange,dataMockup); colormap hot; axis image
    scatter(moleLocs(:,1),moleLocs(:,2),'g');
    for j =1:numMoles
    text(moleLocs(j,1)+10,moleLocs(j,2),num2str(j),'color','w','fontSize',18);
    end
    xlabel('x position (pix)'); ylabel('y position (pix)');
    
	dlg_title = 'What fiducials would you like to average for correction?';
    prompt = {'Enter fid numbers',...
              'What color channel is this in?',...
              };
    def = {num2str(1:numMoles),channel}; % e.g. if 320x320 of 512x512 in far corner, use 193, 193
    num_lines = 1;
    inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
    fidsToUse = eval(['[' inputdialog{1} ']']);
    fidChannel = inputdialog{2};
    close(hSelFig);
    
    %% allow user to use a different spectrum for fid, or use spatially dependent calibration if selected above
    [fidSaveFile, fidSavePath] = uigetfile({'*.mat';'*.*'},...
        'If drift fiducial has different spectrum, open its calibration.mat or easy-dhpsf save now (optional: hit cancel to skip)',...
        'MultiSelect', 'off');
    
    if ~isequal(fidSaveFile,0)
        if ~isequal(fidSaveFile,'calibration.mat')
            load([fidSavePath,fidSaveFile],'s');
            if ~strcmp(s.channel,'0')
                clear s
                calChan = channel;
                load([fidSavePath,fidSaveFile],channel);
            else
                calChan = 's';
            end
            fidCalFile = [eval([calChan '.calFilePrefix']) 'calibration.mat'];
            clear(calChan);
        else
            fidCalFile = [fidSavePath fidSaveFile];
        end
        PSFfits = makeLocalCals(PSFfits,fidCalFile,'FID',0,undoZangleZero);
    elseif spatialCorr
        PSFfits = makeLocalCals(PSFfits,calFile,'FID',0,undoZangleZero);
    end
    
    
    fidTrackX = NaN(size(PSFfits,1)/numMoles,numMoles);
    fidTrackY = NaN(size(PSFfits,1)/numMoles,numMoles);
    fidTrackZ = NaN(size(PSFfits,1)/numMoles,numMoles);
    
    

    for molecule = 1:numMoles
        % extract fitting parameters for this molecule
        moleculeFitParam = PSFfits(PSFfits(:,2) == molecule, :);
        
        % only use good fits as defined by fitting function
        goodFitFlag(:,molecule) = moleculeFitParam(:,13);
        goodFit = goodFitFlag(:,molecule) > 0;
        
        % only use fits that have correct laser on (from sif log)
        if ischar(logFile)
            logFile = cellstr(logFile);
        end
        if ~isequal(logPath,0) % assuming only one sif for whole fiducial track
            % NOTE that this is currently not in use
            sifLogData =  importdata([logPath logFile{1}]);
            %sifLogData = sifLogData(absFrameNum:absFrameNum+numFrames-1,:);
            %absFrameNum = numFrames;
            if channel == 'g'
                goodFit = sifLogData(:,2) == 1 & goodFit;
            elseif channel == 'r'
                goodFit = sifLogData(:,3) == 1 & goodFit;
            end
        end
        
        % raw positions of the fiducial tracks
        fidTrackX(goodFit,molecule) = moleculeFitParam(goodFit,21);
        fidTrackY(goodFit,molecule) = moleculeFitParam(goodFit,22);
        fidTrackZ(goodFit,molecule) = moleculeFitParam(goodFit,23);
        %     numPhotons(:,molecule) = moleculeFitParam(:,17);
        
        
    end
else
    fidTrackX = NaN;
    fidTrackY = NaN;
    fidTrackZ = NaN;
    fidsToUse = [];
    fidChannel = [];
end

%% concatenate data with fid corrections / just concatenate data
if useFidCorrections
    %% compute movement of fiduciaries
    devX = zeros(numFrames,numMoles);
    devY = zeros(numFrames,numMoles);
    devZ = zeros(numFrames,numMoles);
    goodFitFlag = zeros(numFrames,numMoles);
    numPhotons = zeros(numFrames,numMoles);
    avgDevX = zeros(numFrames,1);
    avgDevY = zeros(numFrames,1);
    avgDevZ = zeros(numFrames,1);
    numValidFits = zeros(numFrames,1);
    
    %     textHeader = {'frame number' 'deviation in x (nm)' 'deviation in y (nm)' ...
    %         'deviation in z (nm)' 'good fit flag' 'number of photons'};
    
    syncFrames = zeros(1,numSyncFrames);
    lastGoodFrame = numFrames;
    for a = 1:numSyncFrames % complicated way to find last 25 good frames, to set reference point
        while sum(PSFfits(PSFfits(PSFfits(:,1)==lastGoodFrame,13)>0,2)) ~= numMoles/2*(1+numMoles)
            lastGoodFrame = lastGoodFrame - 1;
        end
        syncFrames(a)=lastGoodFrame;
        lastGoodFrame = lastGoodFrame - 1;
    end
    
    for molecule = fidsToUse
        %% extract fitting parameters for this molecule
        moleculeFitParam = PSFfits(PSFfits(:,2) == molecule, :);
        
        goodFitFlag(:,molecule) = moleculeFitParam(:,13);
        goodFit = goodFitFlag(:,molecule) > 0;
        
        % only use fits that have correct laser on (from sif log)
        % same format as for fidTrackX/Y/Z, but redundant in case that is
        % changed
        if ischar(logFile)
            logFile = cellstr(logFile);
        end
        if ~isequal(logPath,0) % assuming only one sif for whole fiducial track
            sifLogData =  importdata([logPath logFile{1}]);
            %sifLogData = sifLogData(absFrameNum:absFrameNum+numFrames-1,:);
            %absFrameNum = numFrames;
            if channel == 'g'
                goodFit = sifLogData(:,2) == 1 & goodFit;
            elseif channel == 'r'
                goodFit = intersect(sifLogData(:,3) == 1,goodFit);
            end
        end
        
        % compute deviation with respect to bead location averaged over last
        % numSyncFrames frames of the movie
        devX(:,molecule) = moleculeFitParam(:,21) -nanmean(moleculeFitParam(goodFit,21));
        devY(:,molecule) = moleculeFitParam(:,22) -nanmean(moleculeFitParam(goodFit,22));
        devZ(:,molecule) = moleculeFitParam(:,23) -nanmean(moleculeFitParam(goodFit,23));
%         devX(:,molecule) = moleculeFitParam(:,21) ...
%             - mean(moleculeFitParam(any(bsxfun(@eq,moleculeFitParam(:,1), syncFrames),2),21));
%         devY(:,molecule) = moleculeFitParam(:,22) ...
%             - mean(moleculeFitParam(any(bsxfun(@eq,moleculeFitParam(:,1), syncFrames),2),22));
%         devZ(:,molecule) = moleculeFitParam(:,23) ...
%             - mean(moleculeFitParam(any(bsxfun(@eq,moleculeFitParam(:,1), syncFrames),2),23));
%         numPhotons(:,molecule) = moleculeFitParam(:,17);
        
        % write fiduciary data to Excel spreadsheet
        %         xlswrite([saveFilePrefix 'fiduciary deviations.xlsx'], ...
        %             [textHeader; num2cell([(1:numFrames)' devX(:,molecule) devY(:,molecule) ...
        %             devZ(:,molecule) goodFitFlag(:,molecule) numPhotons(:,molecule)])], ...
        %             ['fiduciary ' num2str(molecule)]);
        
        % if particle was fit successfully, add its movement to the average
        avgDevX = avgDevX + goodFit.*devX(:,molecule); % 0s when bad fit
        avgDevY = avgDevY + goodFit.*devY(:,molecule);
        avgDevZ = avgDevZ + goodFit.*devZ(:,molecule);
        numValidFits = numValidFits + goodFit;
    end
    avgDevX = avgDevX./numValidFits; % no fits => /0 => nan
    avgDevY = avgDevY./numValidFits;
    avgDevZ = avgDevZ./numValidFits;
    
end


if useFidCorrections
    %     textHeader = {'frame number' 'deviation in x (nm)' 'deviation in y (nm)' ...
    %         'deviation in z (nm)' 'number of valid fiduciaries'};
    %
    %     xlswrite([saveFilePrefix 'fiduciary deviations.xlsx'], ...
    %         [textHeader; num2cell([(1:numFrames)' tempAvgDevX tempAvgDevY  ...
    %         tempAvgDevZ numValidFits])], 'average of all fiduciaries');
    %
    %     xlswrite([saveFilePrefix 'fiduciary deviations.xlsx'], ...
    %         [{'fiduciary files' 'Number of frames'}; ...
    %         fidFiles num2cell(numFramesInFiles')], ...
    %         'concatenation info');
    
    % load localization data
    %     moleFiles = {};
    
    
    for c = 1:fileNum
        % load data
        load([fitFilePrefix{c} 'molecule fits.mat'],'totalPSFfits','numFrames','calBeadIdx');
        
        if c == 1
            tempPSFfits = totalPSFfits(:,1:27);
        else
            totalPSFfits(:,1) = totalPSFfits(:,1) + sum(numFramesInFiles(1:c-1));
            tempPSFfits = [tempPSFfits; totalPSFfits(:,1:27)];
        end
    end
    totalPSFfits = shiftData(tempPSFfits,'SMACM',imgShiftX,imgShiftY,nmPerPixel);
    
    %% possibly switch calibration fiducial for fits
    if currFidIdx~=calBeadIdx && ~spatialCorr
    dlg_title = 'Use currently selected fiducuial?';
                prompt = { 'Do you want to switch the fiducial to the currently selected one?' };
                def =       { 'No'  };
                questiondialog = questdlg(prompt,dlg_title, def);
                % Handle response
                switch questiondialog
                    case 'Yes'
                        useCurrent=true;
                    case 'No'
                        useCurrent=false;
                    case 'Cancel'
                        error('User cancelled the program');
                end
    else
        useCurrent = false;
    end
    
    if spatialCorr
        [altCalFile, altCalPath] = uigetfile({'*.mat';'*.*'},...
        'If you want to apply another calibration, click on the "calibration.mat" now (optional: hit cancel to skip)',...
        'MultiSelect', 'off');
        if ~isequal(altCalFile,0)
            calFile = [altCalPath altCalFile];
        end
        totalPSFfits = makeLocalCals(totalPSFfits,calFile,'SMACM',0,undoZangleZero);
    elseif useCurrent
        totalPSFfits = makeLocalCals(totalPSFfits,calFile,'SMACM',currFidIdx,undoZangleZero);
        disp('Note that "use current calibration" functionality has not been implemented yet for fiducials');
    end
    
    %     numFrames = sum(numFramesInFiles);
    clear tempPSFfits;
    
    %% De-noise the fiduciary tracks
    
    [avgDevX_denoised,avgDevY_denoised,avgDevZ_denoised] = f_waveletFidTracks(avgDevX,avgDevY,avgDevZ,0);
    
    % apply fiduciary corrections
    if useDenoising == 1
        xFidCorrected = totalPSFfits(:,25) - avgDevX_denoised(totalPSFfits(:,1));
        yFidCorrected = totalPSFfits(:,26) - avgDevY_denoised(totalPSFfits(:,1));
        zFidCorrected = totalPSFfits(:,27) - avgDevZ_denoised(totalPSFfits(:,1));
    else
        xFidCorrected = totalPSFfits(:,25) - avgDevX(totalPSFfits(:,1));
        yFidCorrected = totalPSFfits(:,26) - avgDevY(totalPSFfits(:,1));
        zFidCorrected = totalPSFfits(:,27) - avgDevZ(totalPSFfits(:,1));
    end
    totalPSFfits = [totalPSFfits xFidCorrected yFidCorrected zFidCorrected];
    clear tempAvgDevX tempAvgDevY tempAvgDevZ;
    
    % output corrected data
    %     save([saveFilePrefix 'molecule fits.mat']);
    
    % output excel spreadsheet
    %     textHeader = [textHeader {'fiduciary corrected x location (nm)' ...
    %         'fiduciary corrected y location (nm)' ...
    %         'fiduciary corrected z location (nm)'}];
    
    %     xlswrite([saveFilePrefix 'molecule fits.xlsx'], [textHeader; ...
    %         num2cell(totalPSFfits)], ...
    %         'PSF fits');
    %     xlswrite([saveFilePrefix 'molecule fits.xlsx'], ...
    %         [{'fiduciary files' 'localization files' 'Number of frames'}; ...
    %         fidFiles moleFiles num2cell(numFramesInFiles')], ...
    %         'concatenation info');
    
else
    % just concatenate PSF fits together without compensating for drift
    for fileNum=1:length(fitFilePrefix)
        %         moleFiles = [moleFiles; {[molePath moleFile]}];
        clear numFrames;
        load([fitFilePrefix{fileNum} 'molecule fits.mat'],'totalPSFfits','numFrames','calBeadIdx');
        
        if fileNum == 1
            tempPSFfits = totalPSFfits(:,1:27);
            if ~exist('numFrames','var')
                numFrames = max(totalPSFfits(:,1));
            end
            numFramesInFiles = numFrames;
        else
            if ~exist('numFrames','var')
                numFrames = max(totalPSFfits(:,1));
            end
            numFramesInFiles = [numFramesInFiles numFrames];
            totalPSFfits(:,1) = totalPSFfits(:,1) + sum(numFramesInFiles(1:fileNum-1));
            tempPSFfits = [tempPSFfits; totalPSFfits(:,1:27)];
        end
        
        %         fileNum = fileNum+1;
        %         [moleFile molePath] = uigetfile({'*.mat';'*.*'},...
        %             ['Open data file #' num2str(fileNum) ' with PSF localizations']);
    end
    % call 'numFrames' to match the fiducial-corrected terminology
    numFrames = numFramesInFiles;
    % includes NaNs for fiducial-corrected position fields
    totalPSFfits = [tempPSFfits nan(size(tempPSFfits,1),3)];
    
    %% possibly switch calibration fiducial for fits
    if currFidIdx~=calBeadIdx && ~spatialCorr
    dlg_title = 'Use currently selected fiducuial?';
                prompt = { 'Do you want to switch the calibration fiducial to the currently selected one?' };
                def =       { 'No'  };
                questiondialog = questdlg(prompt,dlg_title, def);
                % Handle response
                switch questiondialog
                    case 'Yes'
                        useCurrent=true;
                    case 'No'
                        useCurrent=false;
                    case 'Cancel'
                        error('User cancelled the program');
                end
    else
        useCurrent = false;
    end
    
    
    
    if spatialCorr
        [altCalFile, altCalPath] = uigetfile({'*.mat';'*.*'},...
        'If you want to apply another calibration, click on the "calibration.mat" now (optional: hit cancel to skip)',...
        'MultiSelect', 'off');
        if ~isequal(altCalFile,0)
            calFile = [altCalPath altCalFile];
        end
        totalPSFfits = makeLocalCals(totalPSFfits,calFile,'SMACM',0,undoZangleZero);
    elseif useCurrent
        totalPSFfits = makeLocalCals(totalPSFfits,calFile,'SMACM',currFidIdx,undoZangleZero);
        disp('Note that "use current calibration" functionality has not been implemented yet for fiducials');
    end
    
    %     numFrames = sum(numFramesInFiles);
    clear tempPSFfits;
    
    fidTrackX = [];
    fidTrackY = [];
    fidTrackZ = [];
    
    % output concatenated data
    %     save([saveFilePrefix 'molecule fits.mat']);
    
    % output excel spreadsheet
    %     xlswrite([saveFilePrefix 'molecule fits.xlsx'], [textHeader; ...
    %         num2cell(totalPSFfits)], ...
    %         'PSF fits');
    %     xlswrite([saveFilePrefix 'molecule fits.xlsx'], ...
    %         [{'localization files' 'Number of frames'}; ...
    %         moleFiles num2cell(numFramesInFiles')], ...
    %         'concatenation info');
    
end

end % end main function

function [PSFfits] = shiftData(PSFfits,type,imgShiftX,imgShiftY,nmPerPixel)
    if strcmp(type,'SMACM')
        xCol = 18;
        yCol = 19;
%         angCol = 20;
%         xOut = 25;
%         yOut = 26;
%         zOut = 27;
%         goodFitCol = 13;
    elseif strcmp(type,'FID')
        xCol = 14;
        yCol = 15;
%         angCol = 16;
%         xOut = 21;
%         yOut = 22;
%         zOut = 23;
%         goodFitCol = 13;
    end
    PSFfits(:,xCol) = PSFfits(:,xCol) + imgShiftX*nmPerPixel;
    PSFfits(:,yCol) = PSFfits(:,yCol) + imgShiftY*nmPerPixel;
end
function [PSFfits] = makeLocalCals(PSFfits,calFile,type,calBeadIdx,undoZangleZero,interpType)
%     interpType = 'linear';
    interpType = 'NN';
    load(calFile,'absLocs','goodFit_f','meanAngles','meanX','meanY','z','zAngleZero');
    numCals = size(absLocs,1);
    
    if strcmp(type,'SMACM')
        xCol = 18;
        yCol = 19;
        angCol = 20;
        xOut = 25;
        yOut = 26;
        zOut = 27;
        goodFitCol = 13;
    elseif strcmp(type,'FID')
        xCol = 14;
        yCol = 15;
        angCol = 16;
        xOut = 21;
        yOut = 22;
        zOut = 23;
        goodFitCol = 13;
    end
    
    goodCal = true(numCals,1);
    for i = 1:numCals 
        if sum(goodFit_f(1,i,:)) < 10 %&& sum(goodFit_b(1,i,:)) < 10
            calDists(:,i) = nan;
            goodCal(i) = false;
            disp(['fiducial # ' num2str(i) ' did not have enough good fits!']);
        end
    end    
    
    if numCals >= 3 && strcmp(type,'SMACM')
    goodX = absLocs(goodCal,1); goodY = absLocs(goodCal,2);
    calHull = convhull(goodX,goodY);
    calHull = [goodX(calHull) goodY(calHull)];
    withinHull = inpolygon(PSFfits(:,xCol),PSFfits(:,yCol),calHull(:,1),calHull(:,2));
    
    figure; plot(calHull(:,1),calHull(:,2))
    hold on; plot(PSFfits(withinHull,xCol),PSFfits(withinHull,yCol),'.g')
    hold on; plot(PSFfits(~withinHull,xCol),PSFfits(~withinHull,yCol),'.r')
    title('Well-defined points, that fall within calibrations');
    else
        withinHull = true(size(PSFfits,1),1);
    end
    
    
    if strcmp(interpType,'NN')
    % length(PSFfits) x numCals... is this too big?
    calDists = sqrt(...
                   (repmat(PSFfits(:,xCol),1,numCals)-repmat(absLocs(:,1)',length(PSFfits),1)).^2 +...
                   (repmat(PSFfits(:,yCol),1,numCals)-repmat(absLocs(:,2)',length(PSFfits),1)).^2);
    
    [~,calNN] = min(calDists,[],2);
    
    
    % this should be rewritten to be an interpolation - this is not a good
    % way to do it long-term
    nearnessScores = ones(size(calDists))./(calDists.^2);
    nearnessScores = nearnessScores./repmat(sum(nearnessScores,2),1,size(nearnessScores,2));

        
    goodFits = PSFfits(:,goodFitCol)>0;
    
    if exist('calBeadIdx') && calBeadIdx ~= 0
        calNN(:)=calBeadIdx;
    end
    %goodFit_forward = logical(squeeze(goodFit_f(1,calBeadIdx,:)));
    
    % for fids: this easily tests whether fid drifts between cals
    % figure; hist(calNN(goodFits,:),1:numCals)
    
    for i = 1:numCals % loop over each cal and redo points near that cal
        if ~goodCal(i)
            continue
        end
        goodFit_forward = logical(squeeze(goodFit_f(1,i,:)));
        PSFfits(calNN==i,xOut) = PSFfits(calNN==i,xCol) ...
            - interp1(squeeze(meanAngles(1,i,goodFit_forward)),...
            squeeze(meanX(1,i,goodFit_forward)),PSFfits(calNN==i,angCol),'spline',nan);
        PSFfits(calNN==i,yOut) = PSFfits(calNN==i,yCol) ...
            - interp1(squeeze(meanAngles(1,i,goodFit_forward)),...
            squeeze(meanY(1,i,goodFit_forward)),PSFfits(calNN==i,angCol),'spline',nan);
        PSFfits(calNN==i,zOut) = interp1(squeeze(meanAngles(1,i,goodFit_forward)),...
            squeeze(z(1,i,goodFit_forward)),PSFfits(calNN==i,angCol),'spline',nan);
        if undoZangleZero
            PSFfits(calNN==i,zOut) = PSFfits(calNN==i,zOut) + zAngleZero(i);
        end
    end
    % end nearest-neighbor module
    elseif strcmp(interpType,'linear')
        % re-extract true observed z values; do not make z=0 at theta=0
        z = squeeze(z(1,1,:))+zAngleZero(1);
        % define vectors repeating across each calibration
        calXrpt = repmat(absLocs(:,1),1,length(z));
        calXrpt = calXrpt+squeeze(meanX); % correct for 'proper drift'
        calXrpt = reshape(calXrpt',[],1);
        
        calYrpt = repmat(absLocs(:,2),1,length(z));
        calYrpt = calYrpt+squeeze(meanY); % correct for 'proper drift'
        calYrpt = reshape(calYrpt',[],1);
        
        calZrpt = repmat(z,numCals,1);
        calTrpt = reshape(squeeze(meanAngles)',[],1);
        goodFitRpt = logical(reshape(squeeze(goodFit_f)',[],1));
        F = scatteredInterpolant(calXrpt(goodFitRpt),calYrpt(goodFitRpt),...
                                 calTrpt(goodFitRpt),calZrpt(goodFitRpt),...
                                 'linear','none');        
        PSFfits(:,zOut) = F(PSFfits(:,xOut),PSFfits(:,yOut),PSFfits(:,angCol));
        % no xy interpolation/correction yet!!
    end % end choice of interpolation modules 
    
    
    
    if exist('withinHull')
        PSFfits(~withinHull,[xOut,yOut,zOut]) = nan; % cannot assign these, since outside range of cals
    end

        
        %% if using nearness scores - be careful with this!!
        % Should be rewritten to be a more careful interpolation
%     PSFfits(:,zOut) = 0;
%     for i = 1:numCals
%         goodFit_forward = logical(squeeze(goodFit_f(1,i,:)));
%         PSFfits(:,zOut) = PSFfits(:,zOut) + nearnessScores(:,i).*interp1(squeeze(meanAngles(1,i,goodFit_forward)),...
%             squeeze(z(1,i,goodFit_forward)),PSFfits(:,angCol),'spline');
%         if undoZangleZero
%         PSFfits(:,zOut) = PSFfits(:,zOut) + nearnessScores(:,i)*zAngleZero(i);
%         end
%     end
end
