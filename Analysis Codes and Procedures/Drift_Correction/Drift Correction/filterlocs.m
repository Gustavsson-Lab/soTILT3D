function [frame_fullfilt,xnm_fullfilt,xinds_fullfilt] = filterlocs(frames,x_nm,framewindow,stdmult,dimension,filepathtemp)
xx = [frames x_nm];
median_xnm = median(x_nm);
mean_xnm = mean(x_nm);
max_xnm = max(x_nm);
min_xnm = min(x_nm);
%framewindow = 1500;
for i = 1:ceil(frames(end)/framewindow)
    framestemp2 = [1:framewindow]+framewindow*(i-1);
    if framestemp2(end) > frames(end)
       framestemp = framestemp2(1):frames(end);
    else
       framestemp = framestemp2;
    end
    frameind2 = (frames>=framestemp(1)& frames<=framestemp(end));
    framestemp3 = frames(frameind2);
    framesegment{i} = frames(frameind2);
    xtemp = x_nm(frameind2);
    maxtemp = max(xtemp);
    mintemp = min(xtemp);
    mediantemp = median(xtemp);
    meantemp = mean(xtemp);
    stdtemp = std(xtemp);
    minmeantemp = abs(maxtemp-mintemp);
    maxmeantemp = abs(maxtemp-mediantemp);
    %x_indtemp{i} = (xtemp>=(mintemp+minmeantemp/1.1)&xtemp<=(maxtemp-maxmeantemp/2));
    x_indtemp{i} = (xtemp>=(meantemp-stdmult*stdtemp)&xtemp<=(meantemp+stdmult*stdtemp));
    x_nmfilt{i} = xtemp(x_indtemp{i});
    frame_filt{i} = framestemp3(x_indtemp{i});
    clear framestemp2 frameind2 framestemp xtemp maxtemp mediantemp meantemp mintemp frametemp3
end
xnm_fullfilt = cell2mat(horzcat(x_nmfilt(:)));
frame_fullfilt = cell2mat(horzcat(frame_filt(:)));
xinds_fullfilt = cell2mat(horzcat(x_indtemp(:)));
figure;
subplot 212
scatter(frame_fullfilt,xnm_fullfilt)
title(['Filtered ',dimension])
%ylim([2000 3000])
subplot 211
scatter(frames,x_nm)
title(['full deviation ',dimension])
saveas(gcf,[filepathtemp,'bead filtered deviation_',dimension,'.jpg'])
saveas(gcf,[filepathtemp,'bead filtered deviation_',dimension,'.fig'])
close(gcf)
%ylim([2000 3000])

% figure;
% subplot 311
% scatter(frame_fullfilt,x_nm(xinds_fullfilt))
% title('x[nm] filtered')
% subplot 312
% scatter(frame_fullfilt,y_nm(xinds_fullfilt))
% title('y[nm] filtered')
% subplot 313
% scatter(frame_fullfilt,z_nm(xinds_fullfilt))
% title('z[nm] filtered')

