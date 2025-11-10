%% Chronological Jinja SNR statistical analysis (excluding February)
% Load SNR observations and their variables
load('merged_arcs.mat', 'snr', 'time', 'elevation', 'azimuth');

%% Turn all cells vertically
all_snr = vertcat(snr{:});
all_time_num = vertcat(time{:});

% Sort time and reorder snr accordingly
[all_time_num_sorted, sortIdx] = sort(all_time_num);
all_snr_sorted = all_snr(sortIdx);

% Convert numeric time to datetime (assuming POSIX time)
all_time = datetime(all_time_num_sorted, 'ConvertFrom', 'posixtime');

%% Remove February from the data
isFeb = month(all_time) == 2;
all_time_nofeb = all_time(~isFeb);
all_snr_nofeb = all_snr_sorted(~isFeb);

%% Downsample raw data for plotting (optional)
ds = 1;  
time_ds = all_time_nofeb(1:ds:end);
snr_ds = all_snr_nofeb(1:ds:end);

%% Compute trendline (moving Mean)
window_size = 5000;   
snr_trend = movmean(all_snr_nofeb, window_size);

% Downsample trendline for plotting (optional)
trend_ds = 1;
time_trend = all_time_nofeb(1:trend_ds:end);
snr_trend_ds = snr_trend(1:trend_ds:end);

%% Main plot (raw + trend)
figure('Position',[100 100 1200 800]);

subplot(3,1,[1 2]); 
plot(time_ds, snr_ds, 'Color', [0.6 0.6 0.9], 'LineWidth', 0.5); 
hold on;
plot(time_trend, snr_trend_ds, 'r-', 'LineWidth', 2);
ylabel('SNR (dB)');
title('SNR Observations at Jinja, Uganda');
grid on;
legend('Raw Data', 'Trendline', 'Location', 'best');

% Add a 10-day buffer at the start and end of x-axis
buffer = days(10);
xlim([min(all_time_nofeb)-buffer, max(all_time_nofeb)+buffer]);

%% =======================
%  OVERALL STATISTICS (no Feb)
% =======================
mean_snr     = mean(all_snr_nofeb, 'omitnan');
median_snr   = median(all_snr_nofeb, 'omitnan');
std_snr      = std(all_snr_nofeb, 'omitnan');
min_snr      = min(all_snr_nofeb);
max_snr      = max(all_snr_nofeb);
p5_snr       = prctile(all_snr_nofeb, 5);
p25_snr      = prctile(all_snr_nofeb, 25);
p75_snr      = prctile(all_snr_nofeb, 75);
p95_snr      = prctile(all_snr_nofeb, 95);
n_obs        = numel(all_snr_nofeb);
time_span    = max(all_time_nofeb) - min(all_time_nofeb);

% --- Print overall stats ---
fprintf('\n===== Overall SNR Statistics =====\n');
fprintf('Obs: %d | Time span: %s\n', n_obs, char(time_span));
fprintf('Mean: %.2f | Median: %.2f | SD: %.2f\n', mean_snr, median_snr, std_snr);
fprintf('Min: %.2f | Max: %.2f\n', min_snr, max_snr);
fprintf('5th: %.2f | 25th: %.2f | 75th: %.2f | 95th: %.2f\n', ...
    p5_snr, p25_snr, p75_snr, p95_snr);

%% =======================
%  DAILY STATS (no Feb) -> only for plot
% =======================
[unique_days, ~, idx_day] = unique(dateshift(all_time_nofeb,'start','day'));
daily_Mean    = accumarray(idx_day, all_snr_nofeb, [], @mean);
daily_std    = accumarray(idx_day, all_snr_nofeb, [], @std);

%% =======================
%  WEEKLY STATS (no Feb)
% =======================
[unique_weeks, ~, idx_week] = unique(dateshift(all_time_nofeb,'start','week'));
weekly_Mean    = accumarray(idx_week, all_snr_nofeb, [], @mean);
weekly_std     = accumarray(idx_week, all_snr_nofeb, [], @std);
weekly_min     = accumarray(idx_week, all_snr_nofeb, [], @min);
weekly_max     = accumarray(idx_week, all_snr_nofeb, [], @max);
weekly_median  = accumarray(idx_week, all_snr_nofeb, [], @median);

% Remove February weeks
isFebWeek = month(unique_weeks) == 2;
unique_weeks(isFebWeek) = [];
weekly_Mean(isFebWeek) = [];
weekly_std(isFebWeek) = [];
weekly_min(isFebWeek) = [];
weekly_max(isFebWeek) = [];
weekly_median(isFebWeek) = [];

% Week labels and start dates
weekLabels = arrayfun(@(i) sprintf('Week %02d', i), 1:numel(unique_weeks), 'UniformOutput', false);
startDates = cellstr(datestr(unique_weeks, 'dd-mmm'));

% Print weekly stats
fprintf('\n===== Weekly SNR Statistics =====\n');
fprintf('%-8s | %-10s | %-6s | %-6s | %-6s | %-6s | %-6s\n', ...
    'Week', 'StartDate', 'Mean', 'Median', 'SD', 'Min', 'Max');
fprintf('-----------------------------------------------------------------\n');

for i = 1:numel(unique_weeks)
    fprintf('%-8s | %-10s | %-6.2f | %-6.2f | %-6.2f | %-6.2f | %-6.2f\n', ...
        weekLabels{i}, startDates{i}, ...
        weekly_Mean(i), weekly_median(i), ...
        weekly_std(i), weekly_min(i), weekly_max(i));
end

%% =======================
%  MONTHLY STATS (no Feb)
% =======================
[unique_months, ~, idx_month] = unique(dateshift(all_time_nofeb,'start','month'));
monthly_Mean    = accumarray(idx_month, all_snr_nofeb, [], @mean);
monthly_std    = accumarray(idx_month, all_snr_nofeb, [], @std);
monthly_min    = accumarray(idx_month, all_snr_nofeb, [], @min);
monthly_max    = accumarray(idx_month, all_snr_nofeb, [], @max);
monthly_median = accumarray(idx_month, all_snr_nofeb, [], @median);

fprintf('\n===== Monthly Statistics =====\n');
fprintf('%-10s | %-6s | %-6s | %-6s | %-6s | %-6s\n', ...
    'Month','Mean','Median','SD','Min','Max');
fprintf('-------------------------------------------------------------\n');
for i = 1:numel(unique_months)
    fprintf('%-10s | %-6.2f | %-6.2f | %-6.2f | %-6.2f | %-6.2f\n', ...
        datestr(unique_months(i), 'mmm yyyy'), ...
        monthly_Mean(i), monthly_median(i), ...
        monthly_std(i), monthly_min(i), monthly_max(i));
end

%% =======================
%  PLOT DAILY, WEEKLY & MONTHLY (no Feb) - SINGLE Y AXIS
% =======================
subplot(3,1,3);
hold on;

%% DAILY (blue line)
h1 = plot(unique_days, daily_Mean, 'b-', 'LineWidth', 1);

%% WEEKLY (magenta line with magenta dots, excl. Feb)
h2 = plot(unique_weeks, weekly_Mean, '-o', 'LineWidth', 1.5, 'Color', 'k', ...
    'MarkerSize', 6, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

%% MONTHLY (green line with black dots, excl. Feb)
h3 = plot(unique_months, monthly_Mean, '-o', 'LineWidth', 1.5, 'Color', 'm', ...
    'MarkerSize', 6, 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'm');

%% Labels & formatting
ylabel('Mean SNR (dB)');
xlabel('Time');
title('Daily, Weekly & Monthly Mean SNR');
grid on;

% Format x-axis with month names
ax = gca;
ax.XTick = unique_months;
ax.XTickLabel = cellstr(datestr(unique_months, 'mmm yyyy'));
ax.XTickLabelRotation = 45;

% Legend
legend([h1, h2, h3], {'Daily Mean', 'Weekly Mean', 'Monthly Mean'}, ...
       'Location','best');
   
% Add a 10-day buffer at the start and end of x-axis
buffer = days(10);
xlim([min(all_time_nofeb)-buffer, max(all_time_nofeb)+buffer]);
