%AVG line scans 
%For LaminAC to LAP2
%Used for soTILT3D manuscript, Sept 27, 2023, NS
%% Information about this code: 

%This code first takes all csv files with a certain prefix and imports them
% into the workspace. Each of the csvs corresponds to the x and y columns of
% the Gaussian fit for each target.
%These csvs can be generated using the code ""
%Then, this code averages the Gaussian fits of all imported fits of for
%each target and fits those curves to a Gaussian
%then finally, this code generates two versions of a boxplot: one with
%outliers and one filtered to be wihtout outliers (percentage based, 5-95%
%confidence interval) 

%line scans 


clc;
clear;

%% Import all Lap 2 x Data

% Define the folder path where your CSV files are located
folder_path = 'Z:\ag134\soTILT3D with PAINT\Setup\FINAL NUCLEUS ANALYSIS\Single 500nm Slice for line scans\single slice Line scans\Nahima line scans\ALL_LINE_SCANS\ALL LINE SCANS - ZEROED TO LAP2';

% Define the prefix of your CSV files
file_prefix = 'LAP2_x_LS';

% Get a list of all CSV files with the same prefix in the folder
csv_files = dir(fullfile(folder_path, [file_prefix '*.csv']));

% Loop through each CSV file and add it to the workspace
for i = 1:numel(csv_files)
    % Read the CSV file into a table
    csv_data = readtable(fullfile(csv_files(i).folder, csv_files(i).name));
    
    % Add the table to the workspace with a unique variable name
    var_name = ['LAP2_x' num2str(i)];
    assignin('base', var_name, csv_data);
end

%% Import all Lap 2 y Data

% Define the folder path where your CSV files are located
folder_path = 'Z:\ag134\soTILT3D with PAINT\Setup\FINAL NUCLEUS ANALYSIS\Single 500nm Slice for line scans\single slice Line scans\Nahima line scans\ALL_LINE_SCANS\ALL LINE SCANS - ZEROED TO LAP2';

% Define the prefix of your CSV files
file_prefix = 'LAP2_y_LS';

% Get a list of all CSV files with the same prefix in the folder
csv_files = dir(fullfile(folder_path, [file_prefix '*.csv']));

% Loop through each CSV file and add it to the workspace
for i = 1:numel(csv_files)
    % Read the CSV file into a table
    csv_data = readtable(fullfile(csv_files(i).folder, csv_files(i).name));
    
    % Add the table to the workspace with a unique variable name
    var_name = ['LAP2_y' num2str(i)];
    assignin('base', var_name, csv_data);
end

%% Import all Lamin A/C x Data

% Define the folder path where your CSV files are located
folder_path = 'Z:\ag134\soTILT3D with PAINT\Setup\FINAL NUCLEUS ANALYSIS\Single 500nm Slice for line scans\single slice Line scans\Nahima line scans\ALL_LINE_SCANS\ALL LINE SCANS - ZEROED TO LAP2';

% Define the prefix of your CSV files
file_prefix = 'LAC_x_LS';

% Get a list of all CSV files with the same prefix in the folder
csv_files = dir(fullfile(folder_path, [file_prefix '*.csv']));

% Loop through each CSV file and add it to the workspace
for i = 1:numel(csv_files)
    % Read the CSV file into a table
    csv_data = readtable(fullfile(csv_files(i).folder, csv_files(i).name));
    
    % Add the table to the workspace with a unique variable name
    var_name = ['LAC_x' num2str(i)];
    assignin('base', var_name, csv_data);
end

%% Import all Lamin A/C y Data

% Define the folder path where your CSV files are located
folder_path = 'Z:\ag134\soTILT3D with PAINT\Setup\FINAL NUCLEUS ANALYSIS\Single 500nm Slice for line scans\single slice Line scans\Nahima line scans\ALL_LINE_SCANS\ALL LINE SCANS - ZEROED TO LAP2';

% Define the prefix of your CSV files
file_prefix = 'LAC_y_LS';

% Get a list of all CSV files with the same prefix in the folder
csv_files = dir(fullfile(folder_path, [file_prefix '*.csv']));

% Loop through each CSV file and add it to the workspace
for i = 1:numel(csv_files)
    % Read the CSV file into a table
    csv_data = readtable(fullfile(csv_files(i).folder, csv_files(i).name));
    
    % Add the table to the workspace with a unique variable name
    var_name = ['LAC_y' num2str(i)];
    assignin('base', var_name, csv_data);
end


%% Turn every single LAP2_x into an array
% Define the prefix you're looking for
desiredPrefix = 'LAP2_x';

% Initialize the number of tables to process
numTables = 25;

% Loop through and convert tables to arrays and assign them to workspace
for i = 1:numTables
    % Construct the variable name
    variableName = [desiredPrefix num2str(i)];
    
    % Check if the variable exists in the workspace
    if evalin('base', ['exist(''' variableName ''', ''var'')'])
        % Use evalin to convert the table to an array
        arrayData = table2array(evalin('base', variableName));
        
        % Create a unique variable name (e.g., 'Array_1', 'Array_2', etc.)
        newArrayName = ['LAP2_x_A_' num2str(i)];
        
        % Assign the array to a variable with the unique name in the workspace
        assignin('base', newArrayName, arrayData);
    end
end

%% Turn every single LAP2_y into an array
% Define the prefix you're looking for
desiredPrefix = 'LAP2_y';

% Initialize the number of tables to process
numTables = 25;

% Loop through and convert tables to arrays and assign them to workspace
for i = 1:numTables
    % Construct the variable name
    variableName = [desiredPrefix num2str(i)];
    
    % Check if the variable exists in the workspace
    if evalin('base', ['exist(''' variableName ''', ''var'')'])
        % Use evalin to convert the table to an array
        arrayData = table2array(evalin('base', variableName));
        
        % Create a unique variable name (e.g., 'Array_1', 'Array_2', etc.)
        newArrayName = ['LAP2_y_A_' num2str(i)];
        
        % Assign the array to a variable with the unique name in the workspace
        assignin('base', newArrayName, arrayData);
    end
end

%% %% Turn every single LAC_x into an array
% Define the prefix you're looking for
desiredPrefix = 'LAC_x';

% Initialize the number of tables to process
numTables = 25;

% Loop through and convert tables to arrays and assign them to workspace
for i = 1:numTables
    % Construct the variable name
    variableName = [desiredPrefix num2str(i)];
    
    % Check if the variable exists in the workspace
    if evalin('base', ['exist(''' variableName ''', ''var'')'])
        % Use evalin to convert the table to an array
        arrayData = table2array(evalin('base', variableName));
        
        % Create a unique variable name (e.g., 'Array_1', 'Array_2', etc.)
        newArrayName = ['LAC_x_A_' num2str(i)];
        
        % Assign the array to a variable with the unique name in the workspace
        assignin('base', newArrayName, arrayData);
    end
end

%% %% %% Turn every single LAC_y into an array
% Define the prefix you're looking for
desiredPrefix = 'LAC_y';

% Initialize the number of tables to process
numTables = 25;

% Loop through and convert tables to arrays and assign them to workspace
for i = 1:numTables
    % Construct the variable name
    variableName = [desiredPrefix num2str(i)];
    
    % Check if the variable exists in the workspace
    if evalin('base', ['exist(''' variableName ''', ''var'')'])
        % Use evalin to convert the table to an array
        arrayData = table2array(evalin('base', variableName));
        
        % Create a unique variable name (e.g., 'Array_1', 'Array_2', etc.)
        newArrayName = ['LAC_y_A_' num2str(i)];
        
        % Assign the array to a variable with the unique name in the workspace
        assignin('base', newArrayName, arrayData);
    end
end


%% %% Sum all the Lap2x data points and then divide by number of vectors to get the average LAP2x
% Define the prefix you're looking for
desiredPrefix = 'LAP2_x_A_';

% Initialize the number of arrays to process
numArrays = 25;  % Change this to the actual number of arrays you have

% Initialize a variable to store the sum
LAP2_x_A_sum = zeros(size(eval([desiredPrefix '1'])));  % Initialize with the size of the first array

% Loop through and calculate the sum
for i = 1:numArrays
    % Construct the variable name
    variableName = [desiredPrefix num2str(i)];
    
    % Check if the variable exists in the workspace
    if evalin('base', ['exist(''' variableName ''', ''var'')'])
        % Use evalin to access the array from the workspace
        currentArray = evalin('base', variableName);
        
        % Add the current array to the running sum
        LAP2_x_A_sum = LAP2_x_A_sum + currentArray;
    end
end

% Now 'LAP2_x_A_sum' contains the sum of the arrays with the specified prefix

LAP2_x_A_avg = LAP2_x_A_sum / 25;

%% %% %% Sum all the Lap2y data points and then divide by number of vectors to get the average LAP2Y
% Define the prefix you're looking for
desiredPrefix = 'LAP2_y_A_';

% Initialize the number of arrays to process
numArrays = 25;  % Change this to the actual number of arrays you have

% Initialize a variable to store the sum
LAP2_y_A_sum = zeros(size(eval([desiredPrefix '1'])));  % Initialize with the size of the first array

% Loop through and calculate the sum
for i = 1:numArrays
    % Construct the variable name
    variableName = [desiredPrefix num2str(i)];
    
    % Check if the variable exists in the workspace
    if evalin('base', ['exist(''' variableName ''', ''var'')'])
        % Use evalin to access the array from the workspace
        currentArray = evalin('base', variableName);
        
        % Add the current array to the running sum
        LAP2_y_A_sum = LAP2_y_A_sum + currentArray;
    end
end

% Now 'LAP2_y_A_sum' contains the sum of the arrays with the specified prefix

LAP2_y_A_avg = LAP2_y_A_sum / 25;

%% Sum all the LACx data points and then divide by number of vectors to get the average LACX
% Define the prefix you're looking for
desiredPrefix = 'LAC_x_A_';

% Initialize the number of arrays to process
numArrays = 25;  % Change this to the actual number of arrays you have

% Initialize a variable to store the sum
LAC_x_A_sum = zeros(size(eval([desiredPrefix '1'])));  % Initialize with the size of the first array

% Loop through and calculate the sum
for i = 1:numArrays
    % Construct the variable name
    variableName = [desiredPrefix num2str(i)];
    
    % Check if the variable exists in the workspace
    if evalin('base', ['exist(''' variableName ''', ''var'')'])
        % Use evalin to access the array from the workspace
        currentArray = evalin('base', variableName);
        
        % Add the current array to the running sum
        LAC_x_A_sum = LAC_x_A_sum + currentArray;
    end
end

% Now 'LAC_x_A_sum' contains the sum of the arrays with the specified prefix

LAC_x_A_avg = LAC_x_A_sum / 25;

%% Sum all the LACy data points and then divide by number of vectors to get the average LACY
% Define the prefix you're looking for
desiredPrefix = 'LAC_y_A_';

% Initialize the number of arrays to process
numArrays = 25;  % Change this to the actual number of arrays you have

% Initialize a variable to store the sum
LAC_y_A_sum = zeros(size(eval([desiredPrefix '1'])));  % Initialize with the size of the first array

% Loop through and calculate the sum
for i = 1:numArrays
    % Construct the variable name
    variableName = [desiredPrefix num2str(i)];
    
    % Check if the variable exists in the workspace
    if evalin('base', ['exist(''' variableName ''', ''var'')'])
        % Use evalin to access the array from the workspace
        currentArray = evalin('base', variableName);
        
        % Add the current array to the running sum
        LAC_y_A_sum = LAC_y_A_sum + currentArray;
    end
end

% Now 'LAC_y_A_sum' contains the sum of the arrays with the specified prefix

LAC_y_A_avg = LAC_y_A_sum / 25;


%% %extract values where curves are maxed for distances
disp(LAP2_x_A_avg(LAP2_y_A_avg == max(LAP2_y_A_avg)));
peak_LAP2 = LAP2_x_A_avg(LAP2_y_A_avg == max(LAP2_y_A_avg));
disp(LAC_x_A_avg(LAC_y_A_avg == max(LAC_y_A_avg)));
peak_LAC = LAC_x_A_avg(LAC_y_A_avg == max(LAC_y_A_avg));

%normalize again for plotting
LAP2_y_A_avg_N = normalize(LAP2_y_A_avg, 'range');
LAC_y_A_avg_N = normalize(LAC_y_A_avg, 'range');

%plotting the figure
figure('Name', 'Line scan for 3 targets');
fontsize = 24;
LW = 4;
set(gcf, 'Color', 'white');
hold on;
plot(LAP2_x_A_avg, LAP2_y_A_avg_N,'color', [0.859 0.008 0.655],'linewidth',LW);
plot(LAC_x_A_avg, LAC_y_A_avg_N,'color', [0.859 0.714 0.008],'linewidth',LW);
xlabel('\color{white} Distance (µm)','fontsize',fontsize)
ylabel('\color{white} Normalized Intensity (a.u.)','fontsize',fontsize)
%xlim([0.1 0.5]);
%ylim([0 1]);
set(text, 'Color', 'w');
set(gca,'Fontsize',20, 'Color', 'k');
set(gcf, 'Color', 'k');
set(gca,'Fontsize',20,'Fontweight','bold', 'XColor', 'w', 'YColor', 'w');
box off;

%graph straight lines where those points are
xline(peak_LAP2, "--", 'color', [0.859 0.008 0.655],'linewidth',LW);
xline(peak_LAC, "--", 'color', [0.859 0.714 0.008],'linewidth',LW);
legend('\color{white} LAP2', '\color{white} Lamin A/C', '', '', '', 'fontsize', 18, 'linewidth', LW, 'Location', 'northeast');
legend boxoff;
hold off;

%% Now fit them to Gaussians

%LAP2

% Create a custom Gaussian model function
gaussianModel = fittype(@(a, b, c, x) a * exp(-((x - b).^2) / (2 * c^2)), 'independent', 'x', 'dependent', 'y');


% Fit the data to the Gaussian model
fitResult_LAP2 = fit(LAP2_x_A_avg(:), LAP2_y_A_avg_N(:), gaussianModel, 'StartPoint', [max(LAP2_y_A_avg_N), mean(LAP2_x_A_avg), std(LAP2_x_A_avg)]);

% Display fitted parameters
disp(fitResult_LAP2);

% Plot the data and the fitted curve
figure('Name', 'Gaussian Fit for Lap 2');
plot(LAP2_x_A_avg, LAP2_y_A_avg_N, 'b.');  % Scatter plot of the data points
hold on;
plot(fitResult_LAP2, 'r');  % Fitted curve in red
xlabel('x');
ylabel('y');
title('Gaussian Curve Fitting');
legend('Data', 'Fitted Curve');

%extract x and y vectors from the fit
x_fit_LAP2 = LAP2_x_A_avg;  
y_fit_LAP2 = feval(fitResult_LAP2, x_fit_LAP2);

%LAC 

% Create a custom Gaussian model function
gaussianModel = fittype(@(a, b, c, x) a * exp(-((x - b).^2) / (2 * c^2)), 'independent', 'x', 'dependent', 'y');

% Fit the data to the Gaussian model
fitResult_LAC = fit(LAC_x_A_avg(:), LAC_y_A_avg_N(:), gaussianModel, 'StartPoint', [max(LAC_y_A_avg_N), mean(LAC_x_A_avg), std(LAC_x_A_avg)]);

% Display fitted parameters
disp(fitResult_LAC);

% Plot the data and the fitted curve
figure('Name', 'Gaussian Fit for Lamin AC');
plot(LAC_x_A_avg, LAC_y_A_avg_N, 'b.');  % Scatter plot of the data points
hold on;
plot(fitResult_LAC, 'r');  % Fitted curve in red
xlabel('x');
ylabel('y');
title('Gaussian Curve Fitting');
legend('Data', 'Fitted Curve');

%extract x and y vectors from the fit
x_fit_LAC = LAC_x_A_avg;  
y_fit_LAC = feval(fitResult_LAC, x_fit_LAC);


%% Find peaks for all the datasets for LAP2
% Define the prefixes for x and y datasets
xPrefix = 'LAP2_x_A_';
yPrefix = 'LAP2_y_A_';

% Define the number of datasets
numDatasets = 25;  % Change this to the actual number of datasets you have

% Initialize an array to store the x-coordinates when y is maxed
xMaxed_LAP2 = zeros(1, numDatasets);

% Loop through datasets
for i = 1:numDatasets
    % Construct variable names for x and y datasets
    xVarName = [xPrefix num2str(i)];
    yVarName = [yPrefix num2str(i)];
    
    % Get x and y data from the workspace
    x_data = evalin('base', xVarName);
    y_data = evalin('base', yVarName);
    
    % Find the maximum y value and its index
    [max_y, maxIndex] = max(y_data);
    
    % Find the corresponding x value
    x_maxed = x_data(maxIndex);
    
    % Store the x-coordinate when y is maxed in the 'xMaxed' array
    xMaxed_LAP2(i) = x_maxed;
    
end

% Now 'xMaxed' contains the x-coordinates when y is maxed for all datasets.

%% %% Find peaks for all the datasets for LAC
% Define the prefixes for x and y datasets
xPrefix = 'LAC_x_A_';
yPrefix = 'LAC_y_A_';

% Define the number of datasets
numDatasets = 25;  % Change this to the actual number of datasets you have

% Initialize an array to store the x-coordinates when y is maxed
xMaxed_LAC = zeros(1, numDatasets);

% Loop through datasets
for i = 1:numDatasets
    % Construct variable names for x and y datasets
    xVarName = [xPrefix num2str(i)];
    yVarName = [yPrefix num2str(i)];
    
    % Get x and y data from the workspace
    x_data = evalin('base', xVarName);
    y_data = evalin('base', yVarName);
    
    % Find the maximum y value and its index
    [max_y, maxIndex] = max(y_data);
    
    % Find the corresponding x value
    x_maxed = x_data(maxIndex);
    
    % Store the x-coordinate when y is maxed in the 'xMaxed' array
    xMaxed_LAC(i) = x_maxed;
    
end

% Now 'xMaxed' contains the x-coordinates when y is maxed for all datasets.



%% % Plot Gaussian Fits for the averages
%extract values where curves are maxed for distances
disp(x_fit_LAP2(y_fit_LAP2 == max(y_fit_LAP2)));
peak_LAP2_fit = x_fit_LAP2(y_fit_LAP2 == max(y_fit_LAP2));
disp(x_fit_LAC(y_fit_LAC == max(y_fit_LAC)));
peak_LAC_fit = x_fit_LAC(y_fit_LAC == max(y_fit_LAC));

%normalize again for plotting
y_fit_LAP2_N = normalize(y_fit_LAP2, 'range');
y_fit_LAC_N = normalize(y_fit_LAC, 'range');

%plotting the figure
figure('Name', 'Line scan for 3 targets -- Gaussian fits');
fontsize = 24;
LW = 4;
set(gcf, 'Color', 'white');
hold on;
plot(x_fit_LAP2, y_fit_LAP2_N,'color', [0.859 0.008 0.655],'linewidth',LW);
plot(x_fit_LAC, y_fit_LAC_N,'color', [0.859 0.714 0.008],'linewidth',LW);
xlabel('\color{white} Distance (µm)','fontsize',fontsize)
ylabel('\color{white} Normalized Intensity (a.u.)','fontsize',fontsize)
xlim([-0.4 0.4]);
%ylim([0 1]);
set(text, 'Color', 'w');
set(gca,'Fontsize',20, 'Color', 'k');
set(gcf, 'Color', 'k');
set(gca,'Fontsize',20,'Fontweight','bold', 'XColor', 'w', 'YColor', 'w');
box off;

%graph straight lines where those points are
xline(peak_LAP2_fit, "--", 'color', [0.859 0.008 0.655],'linewidth',LW);
xline(peak_LAC_fit, "--", 'color', [0.859 0.714 0.008],'linewidth',LW);
legend('\color{white} LAP2', '\color{white} Lamin A/C', '', '', 'fontsize', 18, 'linewidth', LW, 'Location', 'northeast');
legend boxoff;
hold off;

%% %make boxplot figure with no outliers

%Box and whisker plots witout outliers
Distances = (rmoutliers(xMaxed_LAC, 'percentiles', [5 95]))*1000;
Distances_A = transpose(Distances);
fig = figure('Name', 'boxplot for distances');
set(fig, 'Color', 'k');
hold on;
h = boxplot(Distances_A, 'Colors', 'w');
set(h, 'LineWidth', 3);
set(gca,'Fontsize',12,'Fontweight','bold', 'XColor', 'w', 'YColor', 'w');

% Customize the appearance of the plot elements
ax = gca; % Get the current axis handle
ax.Color = 'black'; % Set axis background color to black
ax.XColor = 'white'; % Set x-axis color to white
ax.YColor = 'white'; % Set y-axis color to white

% Set the color of outliers to white
outliers = findobj(gca, 'Tag', 'Outliers');
set(outliers, 'MarkerEdgeColor', 'white', 'MarkerFaceColor', 'white');

% Set the line width of the box plots
for i = 1:numel(h)
    if isa(h(i), 'matlab.graphics.chart.primitive.BoxPlot')
        h(i).Box.LineWidth = 2; % Set line width of the box plots
    end
end

title('Box and Whisker Plot');
% Set custom x-axis labels
xLabels = {'Lamin B1 to LAP2'}; % Replace with your labels
xticks(1:numel(xLabels));
xticklabels(xLabels);
ylabel('\color{white} Distance (nm)','fontsize',16);
%ylabel('\color{white} Distance (µm)','fontsize',fontsize)

disp(median(Distances));
disp(median(xMaxed_LAC));
disp(mean(Distances));
disp(mean(xMaxed_LAC));

%% %% %make boxplot figure with no outliers

%Box and whisker plots with outliers
Distances_noA = xMaxed_LAC*1000;
Distances_AA = transpose(Distances_noA);
fig = figure('Name', 'boxplot for distances with outliers');
set(fig, 'Color', 'k');
hold on;
h = boxplot(Distances_AA, 'Colors', 'w');
set(h, 'LineWidth', 3);
set(gca,'Fontsize',12,'Fontweight','bold', 'XColor', 'w', 'YColor', 'w');

% Customize the appearance of the plot elements
ax = gca; % Get the current axis handle
ax.Color = 'black'; % Set axis background color to black
ax.XColor = 'white'; % Set x-axis color to white
ax.YColor = 'white'; % Set y-axis color to white

% Set the color of outliers to white
outliers = findobj(gca, 'Tag', 'Outliers');
set(outliers, 'MarkerEdgeColor', 'white', 'MarkerFaceColor', 'white');

% Set the line width of the box plots
for i = 1:numel(h)
    if isa(h(i), 'matlab.graphics.chart.primitive.BoxPlot')
        h(i).Box.LineWidth = 2; % Set line width of the box plots
    end
end

title('Box and Whisker Plot');
% Set custom x-axis labels
xLabels = {'Lamin B1 to LAP2'}; % Replace with your labels
xticks(1:numel(xLabels));
xticklabels(xLabels);
ylabel('\color{white} Distance (nm)','fontsize',16);
%ylabel('\color{white} Distance (µm)','fontsize',fontsize)

disp(median(Distances));
disp(median(xMaxed_LAP2));
disp(mean(Distances));
disp(mean(xMaxed_LAP2));
