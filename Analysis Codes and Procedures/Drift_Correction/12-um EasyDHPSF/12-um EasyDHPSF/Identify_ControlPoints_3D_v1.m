% Copyright (c)2013-2018, The Board of Trustees of The Leland Stanford
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

function [] = Identify_ControlPoints_3D_v1()
% This function pairs up corresponding localizations in the two channels

% Key outputs:
% 'matched_cp_reflected' and 'matched_cp_transmitted'
%
% [N/A,N/A,frameStart,frameEnd,meanX,meanY,meanZ,...
% stdX,stdY,stdZ,meanPhotons,numMeasurements,N/A]

% Instrument specific parameters
% nmPerPixel = 125.78;
% roiSize = 270;

if exist('p','var') % if left open from previous run
    delete(p)
end

%% Ask user for relevant datafiles
% these are the localizations for control point candidates in the two
% channels
[reflectedFile, reflectedPath] = uigetfile({'*.mat';'*.*'},'Open data file with raw molecule fits in reflected channel');
if isequal(reflectedFile,0)
    error('User cancelled the program');
end
load([reflectedPath reflectedFile]);
totalPSFfits_reflected = [frameNum, xLoc, yLoc, zLoc, sigmaX, sigmaY, sigmaZ, numPhotons meanBkgnd];
ROIreflected = ROI;
clear frameNum xLoc yLoc zLoc sigmaX sigmaY sigmaZ numPhotons meanBkgnd ROI

[transmittedFile, transmittedPath] = uigetfile({'*.mat';'*.*'},'Open data file with raw molecule fits in transmitted channel');
if isequal(transmittedFile,0)
    error('User cancelled the program');
end
load([transmittedPath transmittedFile]);
totalPSFfits_transmitted = [frameNum, xLoc, yLoc, zLoc, sigmaX, sigmaY, sigmaZ, numPhotons meanBkgnd];
ROItransmitted = ROI;
clear frameNum xLoc yLoc zLoc sigmaX sigmaY sigmaZ numPhotons meanBkgnd ROI

roiSize = max([ROIreflected(3:4) ROItransmitted(3:4)]);

% Load in the logfile containing information of which frames correspond to
% stationary or moving z-positions
[logFile, logPath] = uigetfile({'*.dat';'*.mat';'*.*'},'Open sif log file');
if isequal(logFile,0)
    error('User cancelled the program');
end

% Load in a previous tForm file to serve as an initial guess to match control points
% so that no manual picking is needed
[tFormGuessFile, tFormGuessPath] = uigetfile({'*.mat';'*.*'},'Open previous tForm file as initial guess to match control points; click cancel if no tform file');
if ~isequal(tFormGuessFile,0)
    load([tFormGuessPath tFormGuessFile], 'tform');
    tform_initGuess = tform;
    clear tform
end

%% Find control point candidates
disp('Finding CP candidates...')
[ PSFfits_reflected, PSFfits_transmitted, validFrames, maxNumMeasurement,frameAvgStart,matchLimit] = findCPCandidates(...
    logPath, logFile, totalPSFfits_reflected, totalPSFfits_transmitted );
% PSFfits_reflected & _transmitted: [frameNum xLoc yLoc zLoc sX sY sZ
% signal photons, mean background photons/pix, #beads-in-frame]
disp('CP candidates found.')
save('Identify_ControlPoints_3D_workspace.mat');

%% Calculate the average positions of beads at each stationary point

load('Identify_ControlPoints_3D_workspace.mat');

%[i*ones(beads,1), (1:beads)', frameRange(1)*ones(beads,1), frameRange(2)*ones(beads,1),
%meanX(goodLoc), meanY(goodLoc), meanZ(goodLoc),stdX(goodLoc), stdY(goodLoc),
% stdZ(goodLoc), meanPhotons(goodLoc), numMeasurements(goodLoc)]
disp('Calculating average bead positions...')
[ Locs_reflected, Locs_transmitted ] = avgBeadPos(...
    PSFfits_reflected, PSFfits_transmitted, validFrames, maxNumMeasurement,frameAvgStart);
disp('Average bead positions calculated.')

outputFilePrefix = [reflectedPath filesep 'FilteredLocalizations_std_' ...
    num2str(max(Locs_reflected(:,8))) '_' num2str(max(Locs_reflected(:,9))) '_' num2str(max(Locs_reflected(:,10))) filesep];
mkdir(outputFilePrefix);

% show the averaged bead positions
figure
% set(gcf, 'Position', get(0,'Screensize')); % Maximize figure
scatter3(Locs_reflected(:,5), Locs_reflected(:,6), Locs_reflected(:,7),16, 'filled')
xlabel('x');ylabel('y');zlabel('z');
title('Inspect averaged bead positions in reflected channel. Press any key to continue')
pause
saveas(gcf,[outputFilePrefix 'Locs_reflected_3D.fig']);
figure
% set(gcf, 'Position', get(0,'Screensize')); % Maximize figure
scatter3(Locs_transmitted(:,5), Locs_transmitted(:,6), Locs_transmitted(:,7),16, 'filled')
xlabel('x');ylabel('y');zlabel('z');
title('Inspect averaged bead positions in transmitted channel. Press any key to continue')
pause
saveas(gcf,[outputFilePrefix 'Locs_transmitted_3D.fig']);
close all

save('Identify_ControlPoints_3D_workspace.mat');
save([outputFilePrefix 'Identify_ControlPoints_3D_output.mat']);

%% Find control points in three different z slices
load('Identify_ControlPoints_3D_workspace.mat');

matched_cpLocs_reflected = [];
matched_cpLocs_transmitted = [];
zRanges = [ min(Locs_reflected(:,7)),       max(Locs_reflected(:,7)) ;...
    min(Locs_transmitted(:,7)),     max(Locs_transmitted(:,7)) ];
centerZ = (min(zRanges(:,2))+max(zRanges(:,1)))/2;
zStep = 400;
zSliceLimit = 150;

for zSlice = centerZ-zStep:zStep:centerZ+zStep
    
    % for 2D identification of control points limit the tested range to +- zSliceLimit nm
    
    cpLocs_reflected = Locs_reflected(Locs_reflected(:,7) < (zSlice+zSliceLimit) ...
        & Locs_reflected(:,7) > (zSlice-zSliceLimit),:);
    
    cpLocs_transmitted = Locs_transmitted(Locs_transmitted(:,7) < (zSlice+zSliceLimit) ...
        & Locs_transmitted(:,7) > (zSlice-zSliceLimit),:);
    
    stepCPnum = hist(cpLocs_transmitted(:,1),1:max(cpLocs_transmitted(:,1)));
    
    % only take cpSteps with reasonable number of CP (assuming NHA or
    % similar. Will often fail for instances with only a few CP, and thus
    % no matches.) cpSteps = unique(cpLocs_transmitted(:,1) is equivalent
    % to find(stepCPnum>=1).
    
    cpSteps = find(stepCPnum >= 10);
    
    %     cpSteps = unique(cpLocs_transmitted(:,1));
    
    cpPSFfits_transmitted = PSFfits_transmitted(...
        PSFfits_transmitted(:,4) < (zSlice+zSliceLimit) ...
        & PSFfits_transmitted(:,4) > (zSlice-zSliceLimit),:);
    cpFrames = unique(cpPSFfits_transmitted(:,1));
    
    if exist('tform_initGuess', 'var')
        tform = tform_initGuess;
        
    elseif ~exist('cpChannel1_approx', 'var')   % The first control points need to be hand picked
        [ cpChannel1_approx, cpChannel2_approx, selectedFrame] = ...
            handpickCPs( cpFrames,PSFfits_reflected,PSFfits_transmitted,nmPerPixel,roiSize,ROIreflected,ROItransmitted);
        
        % Take hand-selected pairs and find nearest fit in the respective frame
        [cp_channel1, cp_channel2] = ...
            nearestFit(cpChannel1_approx, cpChannel2_approx,selectedFrame, PSFfits_reflected,PSFfits_transmitted);
        
        % Calculate the (preliminary) transform function based on the
        % handpicked control points using the MATLAB built-in
        % transformation module.
        tform = cp2tform(cp_channel1,cp_channel2,'lwm');
    end
    
    % Use this transform to transform all the average x,y locations in the
    % transmitted channel (Channel 2) to their corresponding location in
    % Channel 1.  (They may not be at the same z in channel 1).
    % Assemble the subset of control points for this slice
    [ matched_cpLocs_reflected_temp, matched_cpLocs_transmitted_temp, matchedCP ] = nearestFitCP(...
        cpLocs_reflected, cpLocs_transmitted, cpSteps, tform,matchLimit );
    
    % This function is evaluated to tell the user how well the
    % control points can be registered within each z-Slice.
    [tform, FRE, TRE, FRE_full, TRE_full] = matlab_transformation(...
        matched_cpLocs_reflected_temp(:,5:6), matched_cpLocs_transmitted_temp(:,5:6), 'affine')
    
    matched_cpLocs_reflected_temp = sortrows(matched_cpLocs_reflected_temp,13);
    matched_cpLocs_reflected_temp(:,13) = matched_cpLocs_reflected_temp(:,13) + size(matched_cpLocs_reflected,1)
    matched_cpLocs_reflected = [matched_cpLocs_reflected ; matched_cpLocs_reflected_temp];
    
    matched_cpLocs_transmitted_temp = sortrows(matched_cpLocs_transmitted_temp,13);
    matched_cpLocs_transmitted_temp(:,13) = matched_cpLocs_transmitted_temp(:,13) + size(matched_cpLocs_transmitted,1)
    matched_cpLocs_transmitted = [matched_cpLocs_transmitted ; matched_cpLocs_transmitted_temp];
    
end

clear FRE FRE_full TRE TRE_full cp_channel1 cp_channel2 zSlice zSliceLimit
clear matchedCP matched_cpLocs_reflected_temp matched_cpLocs_transmitted_temp selectedFrame
save('Identify_ControlPoints_3D_workspace.mat');
save([outputFilePrefix 'Identify_ControlPoints_3D_output.mat']);

%% Evaluate a preliminary 3D transformation
load('Identify_ControlPoints_3D_workspace.mat');

% Assemble the full set of control points
% The parameter fed to this function are chosen empirically, based on
% previous results.
% temp=parcluster;
% matlabpool(temp,temp.NumWorkers-1);
p=parpool('local');
clear temp;
[tform, FRE, TRE, FRE_full, TRE_full] = custom_transformation(...
    matched_cpLocs_reflected(:,5:7),matched_cpLocs_transmitted(:,5:7),'lwquadratic',60,'Gaussian',7,1, true);
delete(p);

cpSteps = unique(Locs_transmitted(:,1));
[ matched_cp_reflected, matched_cp_transmitted, matchedCP ] = nearestFitCP_3D(...
    Locs_reflected, Locs_transmitted, cpSteps, tform);


matched_cp_reflected = sortrows(matched_cp_reflected,13);
matched_cp_transmitted = sortrows(matched_cp_transmitted,13);

save([outputFilePrefix 'Identify_ControlPoints_3D_output.mat']);
disp('3D control point identification complete.')
disp('The key outputs are "matched_cp_reflected" and "matched_cp_transmitted."');
end

%% ---------------------------------------------------------------------------------------------
function [ PSFfits_reflected, PSFfits_transmitted, validFrames, maxNumMeasurement,frameAvgStart,matchLimit] = findCPCandidates(...
    logPath, logFile, totalPSFfits_reflected, totalPSFfits_transmitted )
% This function isolates frames/localizations when there was no xyz motion,
% based on the .sif log. Only filters based on frames.
% These are candidates for control point localizations, stored in
% PSFfits_reflected and PSFfits_transmitted.

%% Ask for user input
dlg_title = 'Please Input Parameters';
nl=sprintf('\n');
prompt = {  'How many stationary frames for each position?',...
    'How many frames to wait before starting the averaging?'...
    ['How closely should the initial transform match points (nm)?' nl...
    '(You may want to use values up to ~800 nm if using an old tform)']...
    };
clear nl
def = {    '50', ...
    '6'...
    '60'...
    };
num_lines = 1;
inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
maxNumMeasurement = str2double(inputdialog{1});
frameAvgStart = str2double(inputdialog{2});
matchLimit = str2double(inputdialog{3});

%% Analyze sif log file
% sifLogData =  importdata([logPath logFile]);
% motionFrames = find(sifLogData(:,1)==-1);   % This finds -1 entries in the shutters correspoding to moving frames
% validPeriod = diff(motionFrames)==maxNumMeasurement+1;
% validMotionFrames = motionFrames(validPeriod);
% validFrames = zeros(length(validMotionFrames),1);
% 
% for i = 1:length(validMotionFrames)
%     frame = validMotionFrames(i); % the frame during which movement occurs
%     temp = (frame+frameAvgStart:frame+maxNumMeasurement)';
%     validFrames = [validFrames; temp];
% end
%%%%%%%%%%

sifLogData = importdata([logPath logFile]);
if size(sifLogData,2) > 2
    logType = 1; % siflog format
    sifLogData(:,4) = sifLogData(:,4).*1000; % convert to nm
    motionFrames = find(sifLogData(:,1)==-1);   % This finds -1 entries in the shutters correspoding to moving frames
    validPeriod = diff(motionFrames)==maxNumMeasurement+1;
    validMotionFrames = motionFrames(validPeriod);
validFrames = []; %%%%     validFrames = zeros(length(validMotionFrames),1);
    
    for i = 1:length(validMotionFrames)
        frame = validMotionFrames(i); % the frame during which movement occurs
        temp = (frame+frameAvgStart:frame+maxNumMeasurement)';
        validFrames = [validFrames; temp];
    end
    
elseif strcmp(logFile(end-2:end), 'mat') 
    logType = 3; % z scan done with PI stage/matlab
    load([logPath logFile],'fVec_cp','nFrames_per_step','numCPScans');
    
    fVec_cp = fVec_cp .* -1; % flip sign
    sifLogData = zeros(length(fVec_cp).*nFrames_per_step,2);
    sifLogData(:,1) = 1:length(fVec_cp).*nFrames_per_step; %col1 is frameNum
    
    zVecRep = [];
    for ii=1:length(fVec_cp)
        zVecRep = [zVecRep; repmat(fVec_cp(ii),[nFrames_per_step,1])];
    end
    sifLogData(:,2) = zVecRep-mean(fVec_cp); %fVec_cp is in microns
    sifLogData(:,2) = sifLogData(:,2).*1000; % convert to nm
    
    motionFrames = [];
    temp = 1;
    for ii= 1:max(sifLogData(:,1))./nFrames_per_step
        motionFrames = [motionFrames; temp];
        temp = temp + nFrames_per_step;
    end
    clear temp
    
    validPeriod = diff(motionFrames)==maxNumMeasurement;%+1;
    validMotionFrames = motionFrames(validPeriod);
validFrames = []; %%%%     validFrames = zeros(length(validMotionFrames),1);
    for i = 1:length(validMotionFrames)
        frame = validMotionFrames(i); % the frame during which movement occurs
        temp = (frame+frameAvgStart:frame+maxNumMeasurement)';
        validFrames = [validFrames; temp];
    end
end 

numFrames = size(sifLogData,1);
PSFfits_reflected = [];
PSFfits_transmitted = [];


%% Only keep localization when there was no xyz motion (based on .sif log)
for frame = 1:numFrames
    if logical(sum(frame == validFrames))
        
        temp_reflected = totalPSFfits_reflected(totalPSFfits_reflected(:,1)==frame,:);
        numBeads = (1:size(temp_reflected,1))';
        temp_reflected = [temp_reflected, numBeads];
        PSFfits_reflected = [PSFfits_reflected; temp_reflected];
        
        temp_transmitted = totalPSFfits_transmitted(totalPSFfits_transmitted(:,1)==frame,:);
        numBeads = (1:size(temp_transmitted,1))';
        temp_transmitted = [temp_transmitted, numBeads];
        PSFfits_transmitted = [PSFfits_transmitted; temp_transmitted];
        %         frame
    end
    
end

end

%% ---------------------------------------------------------------------------------------------
function [ Locs_reflected, Locs_transmitted ] = avgBeadPos(...
    PSFfits_reflected, PSFfits_transmitted, validFrames, framesToAverage,frameStart)
% This function calulates the average positions of beads at each stationary point.
% Control point localizations that are too close together ( currently 100 nm) are discarded.
% Control point localizations can be filtered according to localization
% precision and number of measurements of their localization.
% The function is specific to the current format/sequence of the log file (might need to be changed).
% the averaged localizations are called 'Locs_reflected' and
% 'Locs_transmitted.'
% They are #avgd CP x 12 arrays:

transitionFrames = validFrames([1;find(diff(validFrames)>1)+1]); % the frames after validFrames increases by more than 1

% Reflected Channel
Locs_reflected = [];
windowDiffRef = [];

for i = 1:length(transitionFrames)
    
    % frameRange goes from frame after a jump and is as long as (total # frames
    % acquired on the camera-the frame within that window where averaging starts)
    frameRange = [transitionFrames(i), transitionFrames(i)+(framesToAverage-frameStart)];
    
    temp_reflected = PSFfits_reflected(PSFfits_reflected(:,1)>=frameRange(1) ...
        & PSFfits_reflected(:,1)<=frameRange(2),:);
    
    if size(temp_reflected,1) <= 1 % skip if only one frame in window has any beads
        continue
    else
        
        numBeads = max(temp_reflected(:,10)); % max beads of any frame in window
        meanX = zeros(numBeads,1);
        meanY = zeros(numBeads,1);
        meanZ = zeros(numBeads,1);
        stdX = zeros(numBeads,1);
        stdY = zeros(numBeads,1);
        stdZ = zeros(numBeads,1);
        
        diffX = zeros(numBeads,1);
        diffY = zeros(numBeads,1);
        diffZ = zeros(numBeads,1);
        
        meanPhotons = zeros(numBeads,1);
        numMeasurements = zeros(numBeads,1);
        
        % identify the frames that have the maximum number of beads
        maxBeadFrame = temp_reflected(temp_reflected(:,10)==numBeads,1);
        buffer = 100; % XY radius, nm                       % this parameter might need to be adjustable
        % depending on the precision of the measurement
        % and the density of control points
        beads = 0;
        
        
        for j = 1:numBeads
            
            x = temp_reflected(temp_reflected(:,1)==maxBeadFrame(1) & temp_reflected(:,10)==j,2);
            y = temp_reflected(temp_reflected(:,1)==maxBeadFrame(1) & temp_reflected(:,10)==j,3);
            
            if length(x) == 0
                continue
            else
                % determine how many times each bead was measured at this position
                measurements = (temp_reflected(:,2)-repmat(x,size(temp_reflected,1),1)).^2+...
                    (temp_reflected(:,3)-repmat(y,size(temp_reflected,1),1)).^2 < buffer^2;
                %                     & temp_reflected(:,2)>=x-buffer ...
                %                     & temp_reflected(:,3)<=y+buffer ...
                %                     & temp_reflected(:,3)>=y-buffer;
                numMeasurements(j) = sum(measurements);
                % note: why not just look and see if we have multiple
                % measurements of the same frame number?
                if numMeasurements(j) > framesToAverage-frameStart+1  % two overlapping beads (more beads within the spatial area than frames)
                    continue
                end
                
                % determine statistical parameters for this bead
                meanX(j) = mean(temp_reflected(measurements,2));
                meanY(j) = mean(temp_reflected(measurements,3));
                meanZ(j) = mean(temp_reflected(measurements,4));
                
                stdX(j) = std(temp_reflected(measurements,2));
                stdY(j) = std(temp_reflected(measurements,3));
                stdZ(j) = std(temp_reflected(measurements,4));
                
                tempX = temp_reflected(measurements,2);
                tempY = temp_reflected(measurements,3);
                tempZ = temp_reflected(measurements,4);
                diffX(j) = tempX(end) - tempX(1);
                diffY(j) = tempY(end) - tempY(1);
                diffZ(j) = tempZ(end) - tempZ(1);
                
                meanPhotons(j) = mean(temp_reflected(measurements,8));
                
                beads = beads +1;
            end
            
        end
        goodLoc = find(meanX);
        tempArray = [i*ones(beads,1), (1:beads)',...
            frameRange(1)*ones(beads,1), frameRange(2)*ones(beads,1),...
            meanX(goodLoc), meanY(goodLoc), meanZ(goodLoc),...
            stdX(goodLoc), stdY(goodLoc), stdZ(goodLoc),...
            meanPhotons(goodLoc), numMeasurements(goodLoc)...
            diffX(goodLoc) diffY(goodLoc) diffZ(goodLoc)];
        
        %throw away duplicate entries that may occur due to proximity of
        %two beads
        if ~isempty(tempArray)
            tempArray = sortrows(tempArray,5); % meanX
            temp = diff(tempArray(:,5))==0;
            badFit = logical(zeros(size(tempArray,1),1));
            badFit(1:size(temp)) = temp;
            badFit(1+1:1+size(temp)) = temp;
            tempArray = tempArray(~badFit,:);
            tempArray = sortrows(tempArray,2);
            clear temp badFit
            Locs_reflected = cat(1,Locs_reflected, tempArray(:,1:12));
            windowDiffRef = [windowDiffRef;tempArray(:,13:15)];
        end
        
    end
    
end


% Transmitted Channel
Locs_transmitted = [];
windowDiffTrans = [];
for i = 1:length(transitionFrames)
    
    % frameRange goes from frame after a jump and is as long as (total # frames
    % acquired on the camera-the frame within that window where averaging starts)
    frameRange = [transitionFrames(i), transitionFrames(i)+(framesToAverage-frameStart)];
    
    temp_transmitted = PSFfits_transmitted(PSFfits_transmitted(:,1)>=frameRange(1) ...
        & PSFfits_transmitted(:,1)<=frameRange(2),:);
    
    if size(temp_transmitted,1) <= 1
        continue
    else
        
        numBeads = max(temp_transmitted(:,10));
        meanX = zeros(numBeads,1);
        meanY = zeros(numBeads,1);
        meanZ = zeros(numBeads,1);
        stdX = zeros(numBeads,1);
        stdY = zeros(numBeads,1);
        stdZ = zeros(numBeads,1);
        
        diffX = zeros(numBeads,1);
        diffY = zeros(numBeads,1);
        diffZ = zeros(numBeads,1);
        
        meanPhotons = zeros(numBeads,1);
        numMeasurements = zeros(numBeads,1);
        
        % identify the frames that have the maximum number of beads
        maxBeadFrame = temp_transmitted(temp_transmitted(:,10)==numBeads,1);
        beads = 0;
        
        for j = 1:numBeads
            
            x = temp_transmitted(temp_transmitted(:,1)==maxBeadFrame(1) & temp_transmitted(:,10)==j,2);
            y = temp_transmitted(temp_transmitted(:,1)==maxBeadFrame(1) & temp_transmitted(:,10)==j,3);
            
            if length(x) == 0
                continue
            else
                measurements = (temp_transmitted(:,2)-repmat(x,size(temp_transmitted,1),1)).^2+...
                    (temp_transmitted(:,3)-repmat(y,size(temp_transmitted,1),1)).^2 < buffer^2;
                %                 temp_transmitted(:,2)<=x+buffer ...
                %                     & temp_transmitted(:,2)>=x-buffer ...
                %                     & temp_transmitted(:,3)<=y+buffer ...
                %                     & temp_transmitted(:,3)>=y-buffer;
                numMeasurements(j) = sum(measurements);
                if numMeasurements(j) > framesToAverage-frameStart+1  % two overlapping beads,
                    % these parameters are specific to the
                    % current log file format
                    continue
                end
                
                meanX(j) = mean(temp_transmitted(measurements,2));
                meanY(j) = mean(temp_transmitted(measurements,3));
                meanZ(j) = mean(temp_transmitted(measurements,4));
                
                stdX(j) = std(temp_transmitted(measurements,2));
                stdY(j) = std(temp_transmitted(measurements,3));
                stdZ(j) = std(temp_transmitted(measurements,4));
                
                tempX = temp_transmitted(measurements,2);
                tempY = temp_transmitted(measurements,3);
                tempZ = temp_transmitted(measurements,4);
                diffX(j) = tempX(end) - tempX(1);
                diffY(j) = tempY(end) - tempY(1);
                diffZ(j) = tempZ(end) - tempZ(1);
                
                meanPhotons(j) = mean(temp_transmitted(measurements,8));
                beads = beads +1;
            end
            
        end
        goodLoc = find(meanX);
        tempArray = [i*ones(beads,1), (1:beads)',...
            frameRange(1)*ones(beads,1), frameRange(2)*ones(beads,1),...
            meanX(goodLoc), meanY(goodLoc), meanZ(goodLoc),...
            stdX(goodLoc), stdY(goodLoc), stdZ(goodLoc),...
            meanPhotons(goodLoc), numMeasurements(goodLoc)...
            diffX(goodLoc) diffY(goodLoc) diffZ(goodLoc)];
        %throw away duplicate entries that may occur due to proximity of
        %two beads
        if ~isempty(tempArray)
            tempArray = sortrows(tempArray,5);
            temp = diff(tempArray(:,5))==0;
            badFit = logical(zeros(size(tempArray,1),1));
            badFit(1:size(temp)) = temp;
            badFit(1+1:1+size(temp)) = temp;
            tempArray = tempArray(~badFit,:);
            tempArray = sortrows(tempArray,2);
            clear temp badFit
            Locs_transmitted = cat(1,Locs_transmitted, tempArray(:,1:12));
            windowDiffTrans = [windowDiffTrans ; tempArray(:,13:15)];
        end
    end
    
end

%% Generate diagnostic figures to assess equilibration time of objective
% % what z positions give various shifts?
% shiftTPos = windowDiffTrans(:,3) < 0;
% shiftRPos = windowDiffRef(:,3) < 0;
%
% figure;
% subplot(1,2,1)
% hist(Locs_reflected(shiftRPos,7),40);
% xlabel('z position'); ylabel('counts');
% title('z positions giving shifts < 0 nm for ref channel');
%
% subplot(1,2,2)
% hist(Locs_transmitted(shiftTPos,7),40);
% xlabel('z position'); ylabel('counts');
% title('z positions giving shifts < 0 nm for trans channel');
%
% % what is full distribution of shifts? use CP that fill whole avg window
% fullLengthR = Locs_reflected(:,12) >= 15;
% fullLengthT = Locs_transmitted(:,12) >= 15;
%
% figure;
% subplot(1,2,1);
% hist(windowDiffRef(fullLengthR,3),(max(windowDiffRef(fullLengthR,3))-min(windowDiffRef(fullLengthR,3)))/4);
% xlabel('z shift from begin to end frame (nm)');
% ylabel('counts');
% title('distribution of full-window diffs for reflected channel');
% xlim([-100 100]);
%
% subplot(1,2,2);
% hist(windowDiffTrans(fullLengthT,3),(max(windowDiffTrans(fullLengthT,3))-min(windowDiffTrans(fullLengthT,3)))/4);
% xlabel('z shift from begin to end frame (nm)');
% ylabel('counts');
% title('distribution of full-window diffs for trans channel');
% xlim([-100 100]);

%% Ask for user input

scrsz = get(0,'ScreenSize');
h=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');

subplot(2,4,1)
hist(Locs_reflected(:,8),300)
xlim([0 10])
xlabel('Distance (nm)');
ylabel('Frequency');
title('X Localization Precision');
subplot(2,4,2)
hist(Locs_reflected(:,9),300)
xlim([0 10])
xlabel('Distance (nm)');
ylabel('Frequency');
title('Y Localization Precision');
subplot(2,4,3)
hist(Locs_reflected(:,10),2400)
xlim([0 20])
xlabel('Distance (nm)');
ylabel('Frequency');
title('Z Localization Precision');
subplot(2,4,4)
hist(Locs_reflected(:,12), 30)
xlim([0 framesToAverage])
xlabel('Number of Measurements');
ylabel('Frequency');
title('Number of Measurements');

subplot(2,4,5)
hist(Locs_transmitted(:,8),300)
xlim([0 10])
xlabel('Distance (nm)');
ylabel('Frequency');
title('X Localization Precision');
subplot(2,4,6)
hist(Locs_transmitted(:,9),300)
xlim([0 10])
xlabel('Distance (nm)');
ylabel('Frequency');
title('Y Localization Precision');
subplot(2,4,7)
hist(Locs_transmitted(:,10),2400)
xlim([0 20])
xlabel('Distance (nm)');
ylabel('Frequency');
title('Z Localization Precision');
subplot(2,4,8)
hist(Locs_transmitted(:,12), 30)
xlim([0 framesToAverage])
xlabel('Number of Measurements');
ylabel('Frequency');
title('Number of Measurements');

saveas(gcf,['AveragedCPLocs.fig']);
saveas(gcf,['AveragedCPLocs.png']);

% Restricting the range of localization precisions for the control points
% localizations can ensure that only high quality data points are used
% later on, but it also limits the control point density
dlg_title = 'Please Input Parameters';
prompt = {  'Standard deviation X lower bound',...
    'Standard deviation X upper bound',...
    'Standard deviation Y lower bound',...
    'Standard deviation Y upper bound',...
    'Standard deviation Z lower bound',...
    'Standard deviation z upper bound',...
    'Number of Measurements lower bound',...
    'Number of Measurements upper bound'...
    };
def = {    '0', ...
    '4', ...
    '0', ...
    '4', ...
    '0', ...
    '7', ...
    num2str(framesToAverage-frameStart+1-3),...
    num2str(framesToAverage-frameStart+1)...
    };
num_lines = 1;
inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
maxNumMeasurement = str2double(inputdialog{1});


stdXRange = [str2double(inputdialog{1}) str2double(inputdialog{2})];
stdYRange = [str2double(inputdialog{3}) str2double(inputdialog{4})];
stdZRange = [str2double(inputdialog{5}) str2double(inputdialog{6})];
numMeasurementsRange = [str2double(inputdialog{7}) str2double(inputdialog{8})];

% Apply the chosen Filters
Locs_reflected = Locs_reflected(...
    Locs_reflected(:,8)>stdXRange(1) & ...
    Locs_reflected(:,8)<=stdXRange(2) & ...
    Locs_reflected(:,9)>stdYRange(1) & ...
    Locs_reflected(:,9)<=stdYRange(2) & ...
    Locs_reflected(:,10)>stdZRange(1) & ...
    Locs_reflected(:,10)<=stdZRange(2) & ...
    Locs_reflected(:,12)>=numMeasurementsRange(1) & ...
    Locs_reflected(:,12)<=numMeasurementsRange(2) ,...
    :);

Locs_transmitted = Locs_transmitted(...
    Locs_transmitted(:,8)>stdXRange(1) & ...
    Locs_transmitted(:,8)<=stdXRange(2) & ...
    Locs_transmitted(:,9)>stdYRange(1) & ...
    Locs_transmitted(:,9)<=stdYRange(2) & ...
    Locs_transmitted(:,10)>stdZRange(1) & ...
    Locs_transmitted(:,10)<=stdZRange(2) & ...
    Locs_transmitted(:,12)>=numMeasurementsRange(1) & ...
    Locs_transmitted(:,12)<=numMeasurementsRange(2) ,...
    :);


end % end avgBeadPos

%% ------------------------------------------------------------------------

function [ cpChannel1_approx, cpChannel2_approx, selectedFrame ] = handpickCPs(...
    cpFrames,fits_reflected,fits_transmitted,nmPerPixel,roiSize,ROIreflected,ROItransmitted)
% This function shows frames side by side to let the user handpick control points
% The GUI interface and how many control points need to be selected could be optimized
% Also the user should be allowed to terminate the handpicking himself
%%% TODO: need to make sure that fits are NOT used when canceling:
%%% numAssigned seems to keep increasing, which reduces total pairs chosen
%%% before that routine quits if you change your mind about a pair
scrsz = get(0,'ScreenSize');

% Preallocate for vectors
number_handpicked = 15;
cpChannel1_approx = zeros(number_handpicked,2);
cpChannel2_approx = zeros(number_handpicked,2);
framechannels = zeros(number_handpicked,1); numAssigned = 0;  isOKAY = 0;
selectedFrame = zeros(number_handpicked,1);

for i = 1:10:length(cpFrames)
    goodFitFrame = 0;
    grayBox = zeros(1,4);
    close all
    figure('Position',[(scrsz(3)-1280)/2+1 (scrsz(4)-720)/2 1280 720],'color','w',...
        'Toolbar','figure');
    set(gcf,'DefaultTextFontSize',12,'DefaultAxesFontSize',12);
    showFrameAndFits(fits_reflected,fits_transmitted,cpFrames,i,grayBox,roiSize,nmPerPixel,ROIreflected,ROItransmitted);
    %     subplot(1,2,1)
    %     title({'Is this a good frame?';'Hit enter if yes, click if no'})
    %     subplot(1,2,2)
    %     title({['Assigned Control Points so far: ' num2str(numAssigned)]})
    subplot(1,2,1)
    title({'Is this a good frame?';'Press y if yes, n if no';'Press esc if want to cancel program'})
    subplot(1,2,2)
    title({['Assigned Control Points so far: ' num2str(numAssigned)]})
    
    %     goodFitFrame = waitforbuttonpress;
    %     if goodFitFrame == 1
    %         isOKAY = 1;
    %         while isOKAY == 1
    %             [x,y,isOKAY] = clickPairs(numAssigned,isOKAY);
    %             if isOKAY == 0
    %                 break
    %             end
    %             title('Pick some more!')
    %             cpChannel1_approx(numAssigned+1,:) = [x(1),y(1)];
    %             cpChannel2_approx(numAssigned+1,:) = [x(2),y(2)];
    %             framechannels(numAssigned+1,:) = cpFrames(i);
    %             selectedFrame(numAssigned+1) = cpFrames(i);
    %             grayBox(numAssigned+1,:) = [x(1),y(1),x(2),y(2)];
    %             numAssigned = numAssigned + 1;
    %             showFrameAndFits(fits_reflected,fits_transmitted,cpFrames,i,grayBox,roiSize,nmPerPixel);
    %             subplot(1,2,2)
    %             title({['Assigned Control Points so far: ' num2str(numAssigned)]})
    %         end
    %     elseif numAssigned >= number_handpicked
    %         break
    %     else
    %         continue
    %     end
    
    [~,~,button] = ginput(1);
    if button==121 % y pressed
        goodFitFrame=1;
        isOKAY=1;
        while isOKAY == 1
            [x,y,isOKAY,buttonClick] = clickPairs(numAssigned,isOKAY);
            
            if ~isOKAY & (buttonClick==110) % user selected no - get new frame
                break
            elseif isOKAY & (buttonClick==121) % user selected yes, good click
                title('Pick some more!')
                cpChannel1_approx(numAssigned+1,:) = [x(1),y(1)];
                cpChannel2_approx(numAssigned+1,:) = [x(2),y(2)];
                framechannels(numAssigned+1,:) = cpFrames(i);
                selectedFrame(numAssigned+1) = cpFrames(i);
                grayBox(numAssigned+1,:) = [x(1),y(1),x(2),y(2)];
                numAssigned = numAssigned + 1;
                showFrameAndFits(fits_reflected,fits_transmitted,cpFrames,i,grayBox,roiSize,nmPerPixel,ROIreflected, ROItransmitted);
                subplot(1,2,2)
                title({['Assigned Control Points so far: ' num2str(numAssigned)]})
            else
            end
            
        end
    elseif button == 27 % press esc, terminate program
        break
    elseif numAssigned >= number_handpicked % enough points
        break
    else
        continue
    end
    %%%%%%%%%
    disp(numAssigned)
    %%%%%%%%%
end

end % end handPickCPS

%% -------------------------------------------------------------------------

function showFrameAndFits(...
    PSFfits_reflected,PSFfits_transmitted, cpFrames, frame, grayBox, roiSize, nmPerPixel,ROIreflected,ROItransmitted)
% This function plots the control point localizations  frame-by-frame next to each other

if frame<1 || frame>length(cpFrames)
    return;
end

CCDChipSize = 512;
frameCol = 1; nPhotonsCol = 8; xCenterCol = 2; yCenterCol = 3;

numFits_reflected = size(PSFfits_reflected(PSFfits_reflected(:,1)==cpFrames(frame),:),1);
numFits_transmitted = size(PSFfits_transmitted(PSFfits_transmitted(:,1)==cpFrames(frame),:),1);
numFits = max(numFits_reflected, numFits_transmitted);
markerColors = jet(numFits);

% First subplot (Channel1=Reflected)
temp = sortrows(PSFfits_reflected(PSFfits_reflected(:,frameCol)==cpFrames(frame),:),-nPhotonsCol);
subplot(1,2,1)
hold on
for a = 1:numFits_reflected
    scatter(temp(a,xCenterCol),temp(a,yCenterCol),temp(a,nPhotonsCol)/temp(1,nPhotonsCol)*50,...
        'MarkerFaceColor', markerColors(numFits+1-a,:),...
        'MarkerEdgeColor', markerColors(numFits+1-a,:))
end
% xlim([0 roiSize*nmPerPixel])
% ylim([0 roiSize*nmPerPixel])
xlim([ROIreflected(1)-nmPerPixel ROIreflected(1)+ROIreflected(3)+nmPerPixel])
ylim([ROIreflected(2)-nmPerPixel ROIreflected(2)+ROIreflected(4)+nmPerPixel])
axis square ij
xlabel({'Reflected Channel'; ['Frame: ' num2str(cpFrames(frame))];[num2str(numFits_reflected), ' SM localizations']})
plot([0 roiSize*nmPerPixel], [roiSize*nmPerPixel 0])
plot([0 roiSize*nmPerPixel], [0 roiSize*nmPerPixel])


if sum(grayBox)>0
    for k=1:size(grayBox,1)
        plot(grayBox(k,1),grayBox(k,2),'square','LineWidth',2,'MarkerSize',20)
    end
end
hold off

% Second subplot (Channel2=Transmitted)
temp = sortrows(PSFfits_transmitted(PSFfits_transmitted(:,frameCol)==cpFrames(frame),:),-nPhotonsCol);
subplot(1,2,2)
hold on
for a = 1:numFits_transmitted
    scatter(temp(a,xCenterCol),temp(a,yCenterCol),temp(a,nPhotonsCol)/temp(1,nPhotonsCol)*50,...
        'MarkerFaceColor', markerColors(numFits+1-a,:),...
        'MarkerEdgeColor', markerColors(numFits+1-a,:))
end

% xlim([(CCDChipSize-roiSize)*nmPerPixel CCDChipSize*nmPerPixel])
% ylim([(CCDChipSize-roiSize)*nmPerPixel CCDChipSize*nmPerPixel])
xlim([ROItransmitted(1)-nmPerPixel ROItransmitted(1)+ROItransmitted(3)+nmPerPixel])
ylim([ROItransmitted(2)-nmPerPixel ROItransmitted(2)+ROItransmitted(4)+nmPerPixel])
axis square ij
xlabel({'Transmitted Channel'; ['Frame: ' num2str(cpFrames(frame))];[num2str(numFits_transmitted), ' SM localizations']})
plot([(512-roiSize)*nmPerPixel CCDChipSize*nmPerPixel],...
    [CCDChipSize*nmPerPixel (512-roiSize)*nmPerPixel])
plot([(512-roiSize)*nmPerPixel CCDChipSize*nmPerPixel],...
    [(512-roiSize)*nmPerPixel CCDChipSize*nmPerPixel])

if sum(grayBox)>0
    for k=1:size(grayBox,1)
        plot(grayBox(k,3),grayBox(k,4),'square','LineWidth',2,'MarkerSize',20)
    end
end
hold off

end % end showFrameandFits

%% ------------------------------------------------------------------------

function [x, y, isOKAY, buttonClick] = clickPairs(numAssigned, isOKAY)
% This function lets the user click on control point pairs (Left-->Right)

clear x y
subplot(1,2,1)
title({'Pick two molecules';'(Left plot, then Right plot)'})
[x,y] = ginput(2);
subplot(1,2,1)
cplabel1=text(x(1),y(1),['CP',num2str(numAssigned)],'Color',[0,0,0]);
subplot(1,2,2)
cplabel2=text(x(2),y(2),['CP',num2str(numAssigned)],'Color',[0,0,0]);
subplot(1,2,1)
% title({'Does the following look okay?';'Hit enter if OK, click if NO'})
% isOKAY = waitforbuttonpress();
title({'Does the following look okay?';'Press y if yes, n if no and get new frame,'; 'esc if cancel last pick and stay on this frame'})
[~,~,buttonClick] = ginput(1);
isOKAY = (buttonClick==121);% || (buttonClick==27);
if buttonClick == 27
    delete(cplabel1); delete(cplabel2);
end
end

%% ---------------------------------------------------------------------------------------------
function [cp_channel1, cp_channel2] = nearestFit(cp_channel1_approx, cp_channel2_approx, selectedFrame, c1_allfits, c2_allfits)
% This function takes hand-selected pairs and finds the position of the nearest localization in the respective frame

frameCol = 1; xCenterCol = 2; yCenterCol = 3;

% Channel 1
[numMatch,dim] = size(cp_channel1_approx);
cp_channel1 = zeros(numMatch,dim);

for i=1:numMatch
    frameFits1 = c1_allfits((c1_allfits(:,frameCol) == selectedFrame(i)),:);
    totalDifference = sqrt((cp_channel1_approx(i,1)-frameFits1(:,xCenterCol)).^2 + ...
        (cp_channel1_approx(i,2)-frameFits1(:,yCenterCol)).^2);
    [~,matchedCoord] = min(totalDifference);
    cp_channel1(i,:) = [frameFits1(matchedCoord,xCenterCol),frameFits1(matchedCoord,yCenterCol)];
    clear totalDifference frameFits1 matchedCoord
end

% Channel 2
[numMatch,dim] = size(cp_channel2_approx);
cp_channel2 = zeros(numMatch,dim);

for i=1:numMatch
    frameFits2 = c2_allfits((c2_allfits(:,frameCol) == selectedFrame(i)),:);
    totalDifference = sqrt((cp_channel2_approx(i,1)-frameFits2(:,xCenterCol)).^2 + ...
        (cp_channel2_approx(i,2)-frameFits2(:,yCenterCol)).^2);
    [~,matchedCoord] = min(totalDifference);
    cp_channel2(i,:) = [frameFits2(matchedCoord,xCenterCol),frameFits2(matchedCoord,yCenterCol)];
    clear totalDifference frameFits1 matchedCoord
end

end

%% ---------------------------------------------------------------------------------------------
function [ matched_cpLocs_reflected, matched_cpLocs_transmitted, matchedCP ] = nearestFitCP(...
    cpLocs_reflected, cpLocs_transmitted, cpSteps, tform, matchLimit)
% This function identifies the remaining control point pairs based on the
% structure tform (2D transformation)
% If the target localization is within 60 nm of the transformed
% localization the control point pair is kept.

cpLocs_reflected = [cpLocs_reflected, nan(length(cpLocs_reflected),1)];
cpLocs_transmitted = [cpLocs_transmitted, nan(length(cpLocs_transmitted),1)];
matchedCP = 0;
disp('Generating matches with nearestFitCP over {i} cpSteps');

for i = 1:length(cpSteps)
    i
    temp = cpLocs_transmitted(cpLocs_transmitted(:,1)==cpSteps(i),5:6);
    % Use tform to transform all the average x,y locations in the
    % transmitted channel (Channel 2) to their corresponding location in
    % Channel 1.  (They may not be at the same z in channel 1).
    
    if  tform.ndims_in == 3  % if a previous transformation was used
        temp = [temp zeros(size(temp,1),1)];  % to make dimensionality consistent
        trans_cpLocs_transmitted_temp = transformData(temp,tform);
    else
        trans_cpLocs_transmitted_temp = tforminv(tform, temp);
    end
    
    tempX = cpLocs_reflected(cpLocs_reflected(:,1)==cpSteps(i),5);
    tempY = cpLocs_reflected(cpLocs_reflected(:,1)==cpSteps(i),6);
    
    for j = 1:size(trans_cpLocs_transmitted_temp,1)
        
        residual1 = tempX - trans_cpLocs_transmitted_temp(j,1);
        residual2 = tempY - trans_cpLocs_transmitted_temp(j,2); %+ 3600;
        totalDifference = sqrt((residual1).^2 + (residual2).^2);
        [value,matchedCoord] = min(totalDifference);
        
        if value <= matchLimit   % in units of nm, this parameter might need to be made adjustable.
            matchedCP = matchedCP + 1;
            if isnan(cpLocs_reflected(cpLocs_reflected(:,1)==cpSteps(i) & ...
                    cpLocs_reflected(:,5)==tempX(matchedCoord),13))
                cpLocs_reflected(cpLocs_reflected(:,1)==cpSteps(i) & ...
                    cpLocs_reflected(:,5)==tempX(matchedCoord),13)= matchedCP;
                cpLocs_transmitted(cpLocs_transmitted(:,1)==cpSteps(i) & ...
                    cpLocs_transmitted(:,5)==temp(j,1),13)= matchedCP;
            else
                error('Attempting to assign the same control point twice.  Aborting Execution')
            end
            
        end
    end
end

% Align the matched control points between the two channels
matched_cpLocs_reflected = cpLocs_reflected(~isnan(cpLocs_reflected(:,13)),:);
matched_cpLocs_reflected = sortrows(matched_cpLocs_reflected,13);

matched_cpLocs_transmitted = cpLocs_transmitted(~isnan(cpLocs_transmitted(:,13)),:);
matched_cpLocs_transmitted = sortrows(matched_cpLocs_transmitted,13);

if size(matched_cpLocs_reflected,1) ~= size(matched_cpLocs_transmitted,1)
    error('Unequal number of matched control points')
end

end

%% ---------------------------------------------------------------------------------------------
function [ matched_cpLocs_reflected, matched_cpLocs_transmitted, matchedCP ] = nearestFitCP_3D(...
    cpLocs_reflected, cpLocs_transmitted, cpSteps, tform )
% This function identifies the remaining control point pairs in the entire image domain
% based on the structure tform (3D transformation)
% If the target localization is within 50 nm of the transformed
% localization the control point pair is kept.

cpLocs_reflected = [cpLocs_reflected, nan(length(cpLocs_reflected),1)];
cpLocs_transmitted = [cpLocs_transmitted, nan(length(cpLocs_transmitted),1)];
matchedCP = 0;
disp('Generating matches with nearestFitCP over {i} cpSteps');

for i = 1:length(cpSteps)
    cpSteps(i)
    temp = cpLocs_transmitted(cpLocs_transmitted(:,1)==cpSteps(i),5:7);
    trans_cpLocs_transmitted_temp = transformData(temp,tform);
    tempX = cpLocs_reflected(cpLocs_reflected(:,1)==cpSteps(i),5);
    tempY = cpLocs_reflected(cpLocs_reflected(:,1)==cpSteps(i),6);
    tempZ = cpLocs_reflected(cpLocs_reflected(:,1)==cpSteps(i),7);
    
    for j = 1:size(trans_cpLocs_transmitted_temp,1)
        
        residual1 = tempX - trans_cpLocs_transmitted_temp(j,1);
        residual2 = tempY - trans_cpLocs_transmitted_temp(j,2);
        residual3 = tempZ - trans_cpLocs_transmitted_temp(j,3);
        totalDifference = sqrt((residual1).^2 + (residual2).^2 + (residual3).^2);
        [value,matchedCoord] = min(totalDifference);
        
        if value <= 50   % in units of nm
            matchedCP = matchedCP + 1;
            if isnan(cpLocs_reflected(cpLocs_reflected(:,1)==cpSteps(i) & ...
                    cpLocs_reflected(:,5)==tempX(matchedCoord),13))
                
                cpLocs_reflected(cpLocs_reflected(:,1)==cpSteps(i) & ...
                    cpLocs_reflected(:,5)==tempX(matchedCoord),13)= matchedCP;
                cpLocs_transmitted(cpLocs_transmitted(:,1)==cpSteps(i) & ...
                    cpLocs_transmitted(:,5)==temp(j,1),13)= matchedCP;
                
            else
                error('Attempting to assign the same control point twice.  Aborting Execution')
            end
            
        end
    end
end

% Align the matched control points between the two channels
matched_cpLocs_reflected = cpLocs_reflected(~isnan(cpLocs_reflected(:,13)),:);
matched_cpLocs_reflected = sortrows(matched_cpLocs_reflected,13);

matched_cpLocs_transmitted = cpLocs_transmitted(~isnan(cpLocs_transmitted(:,13)),:);
matched_cpLocs_transmitted = sortrows(matched_cpLocs_transmitted,13);

if size(matched_cpLocs_reflected,1) ~= size(matched_cpLocs_transmitted,1)
    error('Unequal number of matched control points')
end

end

%% ---------------------------------------------------------------------------------------------
function [tform, FRE, TRE, FRE_full, TRE_full] = matlab_transformation(...
    cp_channel1, cp_channel2 , tform_mode)
% This function computes a 2D transformation using the given set of control points
% and also computes the associated FRE and TRE values.
% The code is adapted from "Single-Molecule High Resolution Colocalization of Single
% Probes" by L. Stirling Churchman and James A. Spudich in Cold Spring
% Harbor Protocols (2012), doi: 10.1101/pdb.prot067926
% Modified Definitions of FRE and TRE - Andreas Gahlmann, 20110511

% Transform the data using the cp2tform command
tform = cp2tform(cp_channel1,cp_channel2,tform_mode);

% Calculate the metrices to estimate the error associated with this
% Calculate the fiducial registration error (FRE)
trans_cp_channel2 = tforminv(cp_channel2,tform);

% FRE = sqrt(sum(sum((cp_channel1-trans_cp_channel2).^2))/(length(cp_channel1)));
% FRE_full = ((cp_channel1-trans_cp_channel2).^2)/(length(cp_channel1));

FRE_full = sqrt(sum((cp_channel1-trans_cp_channel2).^2,2));
FRE = mean(FRE_full)

% Calculate the target registration error (TRE)
number_cp = length(cp_channel1); % find the number of control points
% Loop through the control points
disp('Finding TRE, looping over {i} control points');
for i=1:number_cp
    i
    remain_cp = [1:i-1 i+1:number_cp]; % take out that control point
    
    % Calculate the transformation without the ith control point
    tform = cp2tform(cp_channel1(remain_cp,:),cp_channel2(remain_cp,:),tform_mode);
    
    % Transform left out control point with above found transformation
    trans_cp_channel2(i,:) = tforminv(cp_channel2(i,:),tform);
    
end

% TRE = sqrt(sum(sum((cp_channel1 - trans_cp_channel2).^2))/(length(cp_channel1)));
% TRE_full = ((cp_channel1 - trans_cp_channel2).^2)/(length(cp_channel1));
TRE_full = sqrt(sum((cp_channel1-trans_cp_channel2).^2,2));
TRE = mean(TRE_full(TRE_full<100000));

% Restore the full transform function again
tform = cp2tform(cp_channel1,cp_channel2,tform_mode);
channel2_trans = tforminv(cp_channel2,tform);

% show the results
figure
distlimit = 30;

subplot(2,2,1)
scatter(cp_channel1(:,1), cp_channel1(:,2))
title({'Reflected Channel';'Channel 1'})
hold on
scatter(channel2_trans(:,1), channel2_trans(:,2), 10, 'filled')
hold off
axis square
subplot(2,2,2)
scatter(cp_channel2(:,1), cp_channel2(:,2))
title({'Transmitted Channel';'Channel 2'})
axis square

subplot(2,2,3)
hist(FRE_full(FRE_full<=distlimit), 20)
title({['Target Registration Error']; [tform_mode ' Transformation']});
xlabel('Distance (nm)');
ylabel('Frequency');
xlim([0 distlimit]);
legend(['Mean = ' num2str(FRE, 3) ' nm']);

subplot(2,2,4)
hist(TRE_full(TRE_full<distlimit),20)
title({['Fiducial Registration Error']; [tform_mode ' Transformation']});
xlabel('Distance (nm)');
ylabel('Frequency');
xlim([0 distlimit]);
legend(['Mean = ' num2str(TRE, 3) ' nm']);

% saveas(gcf,['z0nm_stats_' tform_mode '.fig']);
% saveas(gcf,['z0nm_stats_' tform_mode '.png']);

end
