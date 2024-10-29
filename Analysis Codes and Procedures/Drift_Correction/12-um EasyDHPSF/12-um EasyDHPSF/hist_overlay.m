
% xLoc = dataSets(2).xLoc_driftCorr;yLoc = dataSets(2).yLoc_driftCorr;
% zLoc = vertcat(featData(:).zLoc);

OSfactor = 2; % oversampling factor relative to pixel size
binSize = nmPerPixel/OSfactor;
xEdges = (min(xRange)-OSfactor/2*binSize):binSize:(max(xRange)+OSfactor/2*binSize);
yEdges = (min(yRange)-OSfactor/2*binSize):binSize:(max(yRange)+OSfactor/2*binSize);


xLoc = vertcat(featData(:,1).xLoc);
yLoc = vertcat(featData(:,1).yLoc);
xPx = floor((xLoc-min(xRange))/binSize)+1;
yPx = floor((yLoc-min(yRange))/binSize)+1;
xPx(xPx<1) = nan;
% zHist = zeros(length(yEdges)-1,length(xEdges)-1);
srHist = zeros(length(yEdges)-1,length(xEdges)-1);
pointIdx = sub2ind(size(srHist),yPx,xPx);
pointIdx(isnan(pointIdx)) = [];
%zRange = [min(zLoc(zLoc>=-1000)) max(zLoc(zLoc<=1000))];
for idx=unique(pointIdx)'
    points = idx == pointIdx;
%     zHist(idx) = median(zLoc(points));
    srHist(idx) = length(xLoc(points));
end
% zHist(brightHist == 0) = NaN;

imwrite(uint8(srHist),['hist_eYFP_' num2str(round(nmPerPixel/OSfactor)) ' nm.tif'])
imwrite(whiteLight,'WL.tif')

xLoc = vertcat(featData(:,2).xLoc);
yLoc = vertcat(featData(:,2).yLoc);
xPx = floor((xLoc-min(xRange))/binSize)+1;
yPx = floor((yLoc-min(yRange))/binSize)+1;
xPx(xPx<1) = nan;
% zHist = zeros(length(yEdges)-1,length(xEdges)-1);
srHist = zeros(length(yEdges)-1,length(xEdges)-1);
pointIdx = sub2ind(size(srHist),yPx,xPx);
pointIdx(isnan(pointIdx)) = [];
%zRange = [min(zLoc(zLoc>=-1000)) max(zLoc(zLoc<=1000))];
for idx=unique(pointIdx)'
    points = idx == pointIdx;
%     zHist(idx) = median(zLoc(points));
    srHist(idx) = length(xLoc(points));
end
% zHist(brightHist == 0) = NaN;

imwrite(uint8(srHist),['hist_PAmChy_' num2str(round(nmPerPixel/OSfactor)) ' nm.tif'])