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

% f_averageCalCurves is a function that can average many calibrations
% together (under development!)

clear all
close all

% save the .mat output by this function into the original calibration
% folder, changing the name of the original and calling the new one
% 'calibration.mat'.

% output 1 x #holes x #steps: not elegant, but this way, can change to an
% arbitrary holes / set of templates without changing the downstream code.
zSamples = -1250:25:1250;	% nm

%% ask user for relevant datafiles

[dataFile, dataPath] = uigetfile({'*.mat';'*.*'},...
    'Open calibration MAT file');
if isequal(dataFile,0)
    error('User cancelled the program');
end
load([dataPath dataFile]);

goodFit_forward = logical(squeeze(goodFit_f(1,1,:)));

templateFrames = interp1(squeeze(meanAngles(1,1,goodFit_forward)),...
            1:length(meanAngles(1,1,goodFit_forward)),-60:30:90,'nearest');

goodFit_b = logical(squeeze(goodFit_b));
goodFit_f = logical(squeeze(goodFit_f));
meanAmpRatio = squeeze(meanAmpRatio);
meanAngles = squeeze(meanAngles);
meanInterlobeDistance = squeeze(meanInterlobeDistance);
meanPhotons = squeeze(meanPhotons);
meanX = squeeze(meanX);
meanY = squeeze(meanY);
stdAmpRatio = squeeze(stdAmpRatio);
stdInterlobeDistance = squeeze(stdInterlobeDistance);
stdX = squeeze(stdX);
stdY = squeeze(stdY);
stddevAngles = squeeze(stddevAngles);
stddevPhotons = squeeze(stddevPhotons);
z = squeeze(z);

numFids = size(meanAngles,1);

% restrict all values to within z=[-1um, 1um]
goodFit_f = goodFit_f & z>=min(zSamples) & z<=max(zSamples);
% quick and dirty filter to eliminate outliers
% goodFit_f = goodFit_f & repmat(meanAngles(:,68)>75 & meanAngles(:,68)<82,[1 size(goodFit_f,2)]);
%     & repmat(z(:,68)>-1200 & z(:,68)<-1200,[1 size(goodFit_f,2)]);

% found that first elements in array were highly negative (~-80), not ~80
% as were the second elements (unwrap issue?): throw out
goodFit_f = goodFit_f & repmat([0 0 0 ones(1,size(goodFit_f,2)-3)],numFids,1);

meanAmpRatio = meanAmpRatio(sum(goodFit_f,2)>2,:);
meanAngles = meanAngles(sum(goodFit_f,2)>2,:);
meanInterlobeDistance = meanInterlobeDistance(sum(goodFit_f,2)>2,:);
meanPhotons = meanPhotons(sum(goodFit_f,2)>2,:);
meanX = meanX(sum(goodFit_f,2)>2,:);
meanY = meanY(sum(goodFit_f,2)>2,:);
stdAmpRatio = stdAmpRatio(sum(goodFit_f,2)>2,:);
stdInterlobeDistance = stdInterlobeDistance(sum(goodFit_f,2)>2,:);
stdX = stdX(sum(goodFit_f,2)>2,:);
stdY = stdY(sum(goodFit_f,2)>2,:);
stddevAngles = stddevAngles(sum(goodFit_f,2)>2,:);
stddevPhotons = stddevPhotons(sum(goodFit_f,2)>2,:);
z = z(sum(goodFit_f,2)>2,:);
goodFit_b = goodFit_b(sum(goodFit_f,2)>2,:);
goodFit_f = goodFit_f(sum(goodFit_f,2)>2,:);

%trim zSamples to make sure we don't look at parts of calibration curves
%that were not measured during actual calibration scan
zSamples = zSamples(zSamples>=min(z(:)) & zSamples<=max(z(:)));


%% compute offsets to overlap curves on top of one another
angleCal = zeros(numFids,length(zSamples));
xCal = zeros(numFids,length(zSamples));
yCal = zeros(numFids,length(zSamples));
ampRatioCal = zeros(numFids,length(zSamples));
interlobeDistCal = zeros(numFids,length(zSamples));
photonsCal = zeros(numFids,length(zSamples));

% sample angle, x, y curves at common z positions in order to compare them
for a = 1:numFids
    angleCal(a,:) = interp1(z(a,goodFit_f(a,:)),meanAngles(a,goodFit_f(a,:)),zSamples,'spline');
    xCal(a,:) = interp1(z(a,goodFit_f(a,:)),meanX(a,goodFit_f(a,:)),zSamples,'spline');
    yCal(a,:) = interp1(z(a,goodFit_f(a,:)),meanY(a,goodFit_f(a,:)),zSamples,'spline');
end
meanAngleAvg = mean(angleCal,1);
meanXAvg = mean(xCal,1);
meanYAvg = mean(yCal,1);

offsetAngle = zeros(1,numFids);
offsetX = zeros(1,numFids);
offsetY = zeros(1,numFids);
for a = 1:numFids
    offsetAngle(a) = lsqnonlin(@(x) angleCal(a,:)-meanAngleAvg-x,mean(angleCal(a,:)-meanAngleAvg));
    offsetX(a) = lsqnonlin(@(x) xCal(a,:)-meanXAvg-x,mean(xCal(a,:)-meanXAvg));
    offsetY(a) = lsqnonlin(@(x) yCal(a,:)-meanYAvg-x,mean(yCal(a,:)-meanYAvg));
end

for a = 1:numFids
    angleCal(a,:) = interp1(z(a,goodFit_f(a,:)),meanAngles(a,goodFit_f(a,:))-offsetAngle(a),zSamples,'spline');
    xCal(a,:) = interp1(z(a,goodFit_f(a,:)),meanX(a,goodFit_f(a,:))-offsetX(a),zSamples,'spline');
    yCal(a,:) = interp1(z(a,goodFit_f(a,:)),meanY(a,goodFit_f(a,:))-offsetY(a),zSamples,'spline');
    ampRatioCal(a,:) = interp1(z(a,goodFit_f(a,:)),meanAmpRatio(a,goodFit_f(a,:)),zSamples,'spline');
    interlobeDistCal(a,:) = interp1(z(a,goodFit_f(a,:)),meanInterlobeDistance(a,goodFit_f(a,:)),zSamples,'spline');
    photonsCal(a,:) = interp1(z(a,goodFit_f(a,:)),meanPhotons(a,goodFit_f(a,:)),zSamples,'spline');
end

% prepare "average" calibration curve for output
meanAngles = permute(mean(angleCal,1),[3 1 2]);
meanX = permute(mean(xCal,1),[3 1 2]);
meanY = permute(mean(yCal,1),[3 1 2]);
meanPhotons = permute(mean(photonsCal,1),[3 1 2]);
meanAmpRatio = permute(mean(ampRatioCal,1),[3 1 2]);
meanInterlobeDistance = permute(mean(interlobeDistCal,1),[3 1 2]);
stddevAngles = permute(std(angleCal,0,1),[3 1 2]);
stdX = permute(std(xCal,0,1),[3 1 2]);
stdY = permute(std(yCal,0,1),[3 1 2]);
stddevPhotons = permute(std(photonsCal,0,1),[3 1 2]);
stdAmpRatio = permute(std(ampRatioCal,0,1),[3 1 2]);
stdInterlobeDistance = permute(std(interlobeDistCal,0,1),[3 1 2]);

z = permute(zSamples,[3 1 2]);
goodFit_b = false([1 1 length(zSamples)]);
goodFit_f = true([1 1 length(zSamples)]);
zAngleZero = zeros(1,1);

%% compare bead calibration curves
%lineColors = distinguishable_colors(size(meanAngles,1));

% angle vs z curves
h=figure;
errorbar(squeeze(z),squeeze(meanAngles),squeeze(stddevAngles));
title([num2str(size(angleCal,1)) ' calibration sources']);
xlabel('z (nm)'); ylabel('angle (deg)');
xlim([-1200 1200]);ylim([-90 90]);
print(h,'-dpng',[dataPath dataFile(1:length(dataFile)-3) ' angle vs z.png']);

% x vs z curves
h=figure;
errorbar(squeeze(z),squeeze(meanX),squeeze(stdX));
title([num2str(size(angleCal,1)) ' calibration sources']);
xlabel('z (nm)'); ylabel('x (nm)');
xlim([-1200 1200]);ylim([-50 50]);
print(h,'-dpng',[dataPath dataFile(1:length(dataFile)-3) ' x vs z.png']);

% y vs z curves
h=figure;
errorbar(squeeze(z),squeeze(meanY),squeeze(stdY));
title([num2str(size(angleCal,1)) ' calibration sources']);
xlabel('z (nm)'); ylabel('y (nm)');
xlim([-1200 1200]);ylim([-50 50]);
print(h,'-dpng',[dataPath dataFile(1:length(dataFile)-3) ' y vs z.png']);

% photons vs z curves
h=figure;
errorbar(squeeze(z),squeeze(meanPhotons),squeeze(stddevPhotons));
title([num2str(size(angleCal,1)) ' calibration sources']);
xlabel('z (nm)'); ylabel('photons detected');
xlim([-1200 1200]);
print(h,'-dpng',[dataPath dataFile(1:length(dataFile)-3) ' photons vs z.png']);

% amplitude ratio vs z curves
h=figure;
errorbar(squeeze(z),squeeze(meanAmpRatio),squeeze(stdAmpRatio));
title([num2str(size(angleCal,1)) ' calibration sources']);
xlabel('z (nm)'); ylabel('Amplitude ratio');
xlim([-1200 1200]);
print(h,'-dpng',[dataPath dataFile(1:length(dataFile)-3) ' amp ratio vs z.png']);

% interlobe distance vs z curves
h=figure;
errorbar(squeeze(z),squeeze(meanInterlobeDistance),squeeze(stdInterlobeDistance));
title([num2str(size(angleCal,1)) ' calibration sources']);
xlabel('z (nm)'); ylabel('Interlobe distance');
xlim([-1200 1200]);
print(h,'-dpng',[dataPath dataFile(1:length(dataFile)-3) ' lobe dist vs z.png']);

%% save data

% repmats are used in order to be able to fit in with data down the road:
% can change to any bead's templates without getting confused by not being
% able to find that bead's data within e.g. meanAngles.

goodFit_b = repmat(goodFit_b,[1,numFids,1]);
goodFit_f = repmat(goodFit_f,[1,numFids,1]);
meanAmpRatio = repmat(meanAmpRatio,[1,numFids,1]);
meanAngles = repmat(meanAngles,[1,numFids,1]);
meanInterlobeDistance = repmat(meanInterlobeDistance,[1,numFids,1]);
meanPhotons = repmat(meanPhotons,[1,numFids,1]);
meanX = repmat(meanX,[1,numFids,1]);
meanY = repmat(meanY,[1,numFids,1]);
stdAmpRatio = repmat(stdAmpRatio,[1,numFids,1]);
stdInterlobeDistance = repmat(stdInterlobeDistance,[1,numFids,1]);
stdX = repmat(stdX,[1,numFids,1]);
stdY = repmat(stdY,[1,numFids,1]);
stddevAngles = repmat(stddevAngles,[1,numFids,1]);
stddevPhotons = repmat(stddevPhotons,[1,numFids,1]);
z = repmat(z,[1,numFids,1]);
zAngleZero = repmat(zAngleZero,[1,numFids,1]);

[dataFile, dataPath] = uiputfile({'*.mat';'*.*'},...
    'Save new calibration MAT file');
if isequal(dataFile,0)
    warning('User did not save the calibration MAT file');
end
save([dataPath dataFile],'meanAngles','meanX','meanY','meanAmpRatio',...
    'meanInterlobeDistance','meanPhotons','stddevAngles','stdX','stdY',...
    'stdAmpRatio','stdInterlobeDistance','stddevPhotons','z','goodFit_b',...
    'goodFit_f','zAngleZero','templateFrames','goodFit_forward');
