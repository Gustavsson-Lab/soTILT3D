%3 nuclear target line scan graphs for 500 nm slice of three targets
%Gaussian fits of line scans exported from ThunderSTORM/ImageJ
%Sep 12, 2023 for soTILT3D project

%data
%LB1 fit
LB1_x = (GaussianFitLB1.X_Fit_Gaussian_no_offset);
LB1_y = normalize(GaussianFitLB1.Fit_Gaussian_no_offset, 'range');
%LAP2 fit
LAP2_x = (GaussianFitLAP2.X_Fit_Gaussian_no_offset);
LAP2_y = normalize(GaussianFitLAP2.Fit_Gaussian_no_offset, 'range');
%LAC fit
LAC_x = (GaussianFitLAC.X_Fit_Gaussian_no_offset);
LAC_y = normalize(GaussianFitLAC.Fit_Gaussian_no_offset, 'range');

%extract values where curves are maxed for distances
disp(LB1_x(LB1_y == max(LB1_y)));
peak_LB1 = LB1_x(LB1_y == max(LB1_y));
disp(LAP2_x(LAP2_y == max(LAP2_y)));
peak_LAP2 = LAP2_x(LAP2_y == max(LAP2_y));
disp(LAC_x(LAC_y == max(LAC_y)));
peak_LAC = LAC_x(LAC_y == max(LAC_y));

%plotting the figure
figure('Name', 'Line scan for 3 targets');
fontsize = 24;
LW = 4;
set(gcf, 'Color', 'white');
hold on;
plot(LB1_x-peak_LAP2, LB1_y,'color', [0.094 0.780 0.769],'linewidth',LW);
plot(LAP2_x-peak_LAP2, LAP2_y,'color', [0.859 0.008 0.655],'linewidth',LW);
plot(LAC_x-peak_LAP2, LAC_y,'color', [0.859 0.714 0.008],'linewidth',LW);
xlabel('\color{white} Distance (µm)','fontsize',fontsize)
ylabel('\color{white} Normalized Intensity (a.u.)','fontsize',fontsize)
xlim([-0.3 0.3]);
%ylim([0 1]);
set(text, 'Color', 'w');
set(gca,'Fontsize',20, 'Color', 'k');
set(gcf, 'Color', 'k');
set(gca,'Fontsize',20,'Fontweight','bold', 'XColor', 'w', 'YColor', 'w');
box off;

%graph straight lines where those points are
xline(peak_LB1-peak_LAP2, "--", 'color', [0.094 0.780 0.769],'linewidth',LW);
xline(peak_LAP2-peak_LAP2, "--", 'color', [0.859 0.008 0.655],'linewidth',LW);
xline(peak_LAC-peak_LAP2, "--", 'color', [0.859 0.714 0.008],'linewidth',LW);
legend('\color{white} Lamin B1', '\color{white} LAP2', '\color{white} Lamin A/C', '', '', '', 'fontsize', 18, 'linewidth', LW, 'Location', 'northeast');
legend boxoff;
hold off;

