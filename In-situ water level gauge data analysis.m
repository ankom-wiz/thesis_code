%% =====================================================================
%  GAUGE DATA ANALYSIS TOOL
%  Computes overall, monthly, and weekly statistics for CSV water level
%  gauge data, and also computes high-resolution deviation from daily mean
%  for use in surface roughness or hydrological modeling.
%
%  Additionally, plots the raw water level data vs. time.
%
%  Output:
%   - Command Window summary (overall/monthly/weekly stats)
%   - Plot of in-situ water level gauge
%   - highResGaugeData.mat (per-observation deviation dataset)
%  =====================================================================

%% --- Load data ---
filename = 'gauge_data.csv';
opts = detectImportOptions(filename, 'NumHeaderLines', 0);
opts = setvartype(opts, {'datetime', 'double'}); % 1st col: Date, 2nd col: Gauge reading
data = readtable(filename, opts);

% Name columns for clarity
data.Properties.VariableNames = {'Date', 'Gauge'};

%% --- Basic checks ---
fprintf('\n=== Gauge Data Summary ===\n');
fprintf('Total records: %d\n', height(data));
fprintf('Date range: %s to %s\n', datestr(min(data.Date)), datestr(max(data.Date)));
fprintf('----------------------------------------------------------\n');

%% --- Overall statistics ---
overall.meanVal   = mean(data.Gauge, 'omitnan');
overall.medianVal = median(data.Gauge, 'omitnan');
overall.stdVal    = std(data.Gauge, 'omitnan');
overall.minVal    = min(data.Gauge);
overall.maxVal    = max(data.Gauge);
overall.rangeVal  = overall.maxVal - overall.minVal;

fprintf('\n=== Overall Statistics ===\n');
fprintf('Mean    | Median  | Stdev   | Min     | Max     | Range\n');
fprintf('----------------------------------------------------------\n');
fprintf('%.3f | %.3f | %.3f | %.3f | %.3f | %.3f\n', ...
    overall.meanVal, overall.medianVal, overall.stdVal, ...
    overall.minVal, overall.maxVal, overall.rangeVal);

%% --- Monthly statistics ---
data.Month = month(data.Date);
monthNames = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

monthlyStats = groupsummary(data, 'Month', {'mean','median','std','min','max'}, 'Gauge');
monthlyStats.Range_Gauge = monthlyStats.max_Gauge - monthlyStats.min_Gauge;

fprintf('\n=== Monthly Statistics ===\n');
fprintf('Month | Mean    | Median  | Stdev   | Min     | Max     | Range\n');
fprintf('-----------------------------------------------------------------\n');

for i = 1:height(monthlyStats)
    m = monthlyStats.Month(i);
    fprintf('%s   | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f\n', ...
        monthNames{m}, ...
        monthlyStats.mean_Gauge(i), ...
        monthlyStats.median_Gauge(i), ...
        monthlyStats.std_Gauge(i), ...
        monthlyStats.min_Gauge(i), ...
        monthlyStats.max_Gauge(i), ...
        monthlyStats.Range_Gauge(i));
end

%% --- Weekly statistics ---
[unique_weeks, ~, idx_week] = unique(dateshift(data.Date, 'start', 'week'));

weekly_Mean    = accumarray(idx_week, data.Gauge, [], @mean, NaN);
weekly_std     = accumarray(idx_week, data.Gauge, [], @std, NaN);
weekly_min     = accumarray(idx_week, data.Gauge, [], @min, NaN);
weekly_max     = accumarray(idx_week, data.Gauge, [], @max, NaN);
weekly_median  = accumarray(idx_week, data.Gauge, [], @median, NaN);

weekLabels = arrayfun(@(i) sprintf('Week %02d', i), 1:numel(unique_weeks), 'UniformOutput', false);

fprintf('\n===== Weekly Gauge Statistics =====\n');
fprintf('%-8s | %-10s | %-6s | %-6s | %-6s | %-6s | %-6s\n', ...
    'Week','StartDate','Mean','Median','SD','Min','Max');
fprintf('--------------------------------------------------------------------------\n');

for i = 1:numel(unique_weeks)
    weekStartStr = char(unique_weeks(i), 'dd-MMM');
    fprintf('%-8s | %-10s | %-6.2f | %-6.2f | %-6.2f | %-6.2f | %-6.2f\n', ...
        weekLabels{i}, weekStartStr, ...
        weekly_Mean(i), weekly_median(i), ...
        weekly_std(i), weekly_min(i), weekly_max(i));
end

%% =====================================================================
%  HIGH-RESOLUTION DEVIATION ANALYSIS (from jinja_gauge_highres.m)
% =====================================================================

fprintf('\n=== High-Resolution Daily Deviation Analysis ===\n');

% --- Sort and clean ---
T = sortrows(data, 'Date');
T = rmmissing(T, 'DataVariables', 'Gauge');

% --- Compute daily mean for each observation ---
dayVec = dateshift(T.Date, 'start', 'day');
[G, ~] = findgroups(dayVec);
dailyMean = splitapply(@mean, T.Gauge, G);
dailyMeanObs = dailyMean(G);

% --- Compute deviation per observation ---
deviation = T.Gauge - dailyMeanObs;

% --- Save to MAT file ---
highResGaugeData.Timestamp = T.Date;
highResGaugeData.Deviation = deviation;

save('highResGaugeData.mat', 'highResGaugeData');
fprintf('Saved high-resolution deviation dataset: highResGaugeData.mat\n');

%% =====================================================================
%  PLOT: In-situ Water Level Gauge
% =====================================================================

figure('Color','w');
plot(data.Date, data.Gauge, 'b.-', 'LineWidth', 1.2, 'MarkerSize', 10);
xlabel('Time');
ylabel('Water Level (m)');
title('In-situ Water Level Gauge Measurements');
legend('Water gauge level (m)', 'Location', 'best');
grid on;
datetick('x','mmm','keeplimits'); % Show month names on X-axis
set(gca, 'FontSize', 11);

fprintf('Plot generated: In-situ water level gauge vs. time.\n');
fprintf('===============================================================\n');

%% =====================================================================
%  PLOT: Weekly Mean Water Levels (same X-axis style as monthly plot)
% =====================================================================

figure('Color','w'); % Open in a new window
plot(unique_weeks, weekly_Mean, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Time');
ylabel('Weekly Mean Water Level (m)');
title('Weekly Mean Water Level Gauge Measurements');
grid on;

% --- Match the same month-based X-axis style ---
datetick('x','mmm','keeplimits'); % Label axis with month abbreviations
xlim([min(data.Date) max(data.Date)]); % Match full date range of dataset

set(gca, 'FontSize', 11);
legend('Weekly Mean', 'Location', 'best');
