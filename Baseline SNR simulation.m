%% MPSim baseline simulation (simulate SNR as close as possible to reality)
%% Load GNSS merged arcs (observed data file)
% Contains: in-situ SNR obs. and corresponsing elevation, azimuth, timestamp
load('merged_arcs.mat', 'snr', 'elevation', 'azimuth', 'time');

%% Fix merged_mat parameter data forms
% Concatenate all cells vertically
all_elev = vertcat(elevation{:});
all_azim = vertcat(azimuth{:});
all_time = vertcat(time{:});
all_snr  = vertcat(snr{:});

% Sort by time
[all_time_sorted, sortIdx] = sort(all_time);
all_elev_sorted = all_elev(sortIdx);
all_azim_sorted = all_azim(sortIdx);
all_snr_sorted  = all_snr(sortIdx);

% Convert to datetime
time_dt = datetime(all_time_sorted, 'ConvertFrom', 'posixtime');

% Remove February (only contained 2 days -- incomplete dataset)
isFeb = month(time_dt) == 2;
all_time_nofeb = all_time_sorted(~isFeb);
all_elev_nofeb = all_elev_sorted(~isFeb);
all_azim_nofeb = all_azim_sorted(~isFeb);
all_snr_nofeb  = all_snr_sorted(~isFeb);
time_dt_nofeb  = time_dt(~isFeb);

%% Downsampling option (set to 1 to keep all data)
ds = 1;   % downsampling factor (1 = no downsample)
time_ds = all_time_nofeb(1:ds:end);
elev_ds = all_elev_nofeb(1:ds:end);
azim_ds = all_azim_nofeb(1:ds:end);
snr_ds  = all_snr_nofeb(1:ds:end);
time_dt_ds = time_dt_nofeb(1:ds:end);

%% Water surface statistics from high-resolution gauge data
% This segment assigns the highest possible temporal surface deviation to 
% the model dataset. It segments the daily data into equal parts as many as the
% deviations calculated from gauge data water level temporal variance.
load('highResGaugeData.mat', 'highResGaugeData');  % Has: Timestamp, Deviation

obsDates = dateshift(time_dt_ds, 'start', 'day');    % Convert timestamps: '2024-03-20 00:00:01' -> '2024-03-20'
gaugeDates = dateshift(highResGaugeData.Timestamp, 'start', 'day');  % Gauge dates

numObs = numel(obsDates);
wl_std_vec = nan(numObs,1);
uniqueDays = unique(obsDates);

for i = 1:numel(uniqueDays)
    modelIdx = find(obsDates == uniqueDays(i));   % Indices in model timestamps for this day
    gaugeIdx = find(gaugeDates == uniqueDays(i)); % Indices in high-res deviations for this day
    
    if isempty(gaugeIdx)
        continue;  % Skip if no gauge data for the day
    end
    
    dev_day = highResGaugeData.Deviation(gaugeIdx);
    n_model = numel(modelIdx);
    n_dev = numel(dev_day);
    
    % Determine chunk sizes for each deviation
    chunk_sizes = floor(n_model / n_dev) * ones(n_dev,1);    % base size
    remainder = n_model - sum(chunk_sizes);                % leftover timestamps
    chunk_sizes(1:remainder) = chunk_sizes(1:remainder) + 1; % distribute remainder
    
    % Assign deviations to model timestamps
    idx_start = 1;
    for j = 1:n_dev
        idx_end = idx_start + chunk_sizes(j) - 1;
        wl_std_vec(modelIdx(idx_start:idx_end)) = dev_day(j);
        idx_start = idx_end + 1;
    end
end
%wl_std_vec = abs(wl_std_vec);

%% MPSim settings
%% Step 1: Load default settings
sett = snr_settings();

%% Step 2: Customize antenna
sett.ref.height_ant = 2.6;                  % Antenna height (m)
sett.ant.slope = 'tipped';                  % Antenna tilted 90°
sett.ant.aspect = 176;                      % Azimuth mean = 176°; reflections mainly from southward azimuths
sett.ant.model = 'isotropic horizontal';    % Use patch-like antenna model
sett.ant.radome = 'NONE';                   % Antenna (isotropic model) has no radome

%% Step 3: Define surface properties
% Surface geometry
sett.sfc.fnc_snr_setup_sfc_geometry = @snr_setup_sfc_geometry_horiz; 
sett.sfc.slope = 0;                       % Horizontal
sett.sfc.aspect = 176;                    % Facing south (matches river direction)
sett.sfc.dir_nrml = [0;0;1];              % north-east-up local vector for horizontal surface
sett.sfc.height_std = wl_std_vec;         % Daily water level deviations from gauge data)
sett.sfc.vert_datum = 'top';              % Antenna height is referenced to top of the water
sett.ref.pos_origin = [0, 0, 0];          % Base of antenna pole (local origin)
sett.sfc.pos_sfc0 = sett.ref.pos_origin;  % Intersection of top surface with antenna pole

% Reflecting surface model
% HALFSPACE Surface model
% Simplistic model that only models the signal interaction with air and water
% surface considering the GPS (radio wave) permeability in water
sett.sfc.fnc_snr_setup_sfc_material = @snr_setup_sfc_material_halfspaces;
sett.sfc.material_top = 'air';           % top halfspace = air
sett.sfc.material_bottom = 'freshwater'; % bottom halfspace = water (interface is air-water)

%% Step 4: GNSS and signal configuration
sett.opt.gnss_name = 'GPS';
sett.opt.freq_name = 'L1';
sett.opt.rec.bandwidth_noise = 1;            % Hz
sett.opt.rec.temperature_noise = 300;        % Rec. thermal noise (Kelvin)- room temp. (>290 due to local climate)
sett.opt.rec.ant_density_noise_db = -199.3;  % dB-W/Hz -- thermal environmental noise (empirical testing)

%% Step 5: Observation directions
sett.sat.num_obs = numel(elev_ds);  % Match number of real observations
sett.sat.elev = elev_ds(:);
sett.sat.azim = azim_ds(:);
sett.sat.epoch = time_ds(:);

%% Step 6: Prepare the simulation setup
setup = snr_setup(sett);

%% Step 7: Run simulation (MPSim)
[result, snr_db, phasor_interf, phasor_composite] = snr_fwd(setup);

%% ===============================
%  OVERALL STATISTICS (Simulated only)
% ===============================
snr_data = snr_db;
n_obs = numel(snr_data);
time_span = max(time_dt_ds) - min(time_dt_ds);
mean_snr = mean(snr_data, 'omitnan');
median_snr = median(snr_data, 'omitnan');
std_snr = std(snr_data, 'omitnan');
min_snr = min(snr_data);
max_snr = max(snr_data);
p5_snr = prctile(snr_data,5);
p25_snr = prctile(snr_data,25);
p75_snr = prctile(snr_data,75);
p95_snr = prctile(snr_data,95);

fprintf('\n===== Simulated SNR Statistics =====\n');
fprintf('Obs: %d | Time span: %s\n', n_obs, string(time_span));
fprintf('Mean: %.2f | Median: %.2f | SD: %.2f\n', ...
    mean_snr, median_snr, std_snr);
fprintf('Min: %.2f | Max: %.2f\n', min_snr, max_snr);
fprintf('5th: %.2f | 25th: %.2f | 75th: %.2f | 95th: %.2f\n', ...
    p5_snr, p25_snr, p75_snr, p95_snr);

%% ===============================
%  MONTHLY STATISTICS (Simulated only)
% ===============================
months = dateshift(time_dt_ds,'start','month');
[unique_months, ~, idx_month] = unique(months);

snr_data = snr_db;
monthly_mean   = accumarray(idx_month, snr_data, [], @mean);
monthly_median = accumarray(idx_month, snr_data, [], @median);
monthly_std    = accumarray(idx_month, snr_data, [], @std);
monthly_min    = accumarray(idx_month, snr_data, [], @min);
monthly_max    = accumarray(idx_month, snr_data, [], @max);

fprintf('\n===== Monthly Simulated SNR =====\n');
fprintf('%-10s | %-6s | %-6s | %-6s | %-6s | %-6s\n', ...
    'Month','Mean','Median','SD','Min','Max');
fprintf('---------------------------------------------------\n');
for i = 1:numel(unique_months)
    fprintf('%-10s | %-6.2f | %-6.2f | %-6.2f | %-6.2f | %-6.2f\n', ...
        datestr(unique_months(i),'mmm yyyy'), ...
        monthly_mean(i), monthly_median(i), ...
        monthly_std(i), monthly_min(i), monthly_max(i));
end

%% ===============================
%  DAILY MEAN TRENDLINE PLOTS
% ===============================
days = dateshift(time_dt_ds,'start','day');
[unique_days, ~, idx_days] = unique(days);

daily_mean_obs = accumarray(idx_days, snr_ds, [], @mean);
daily_mean_sim = accumarray(idx_days, snr_db, [], @mean);

figure('Name','Observed vs Simulated SNR (dB)');
plot(unique_days, daily_mean_obs, 'b-', 'LineWidth', 1.5, 'DisplayName','Observed'); hold on;
plot(unique_days, daily_mean_sim, 'r-', 'LineWidth', 1.5, 'DisplayName','Simulated');
xlabel('Date'); ylabel('Daily Mean SNR (dB)');
title('Observed vs Simulated SNR (dB)');
legend('Location','best'); grid on;

%% ===============================
%% Observed vs Simulated SNR (V/V) for individual arcs
%% ===============================

% Plotting sampled individual satellite passes (singular arcs) comparing 
% observed SNR with simulated interferometric voltage as a function of 
% the sine of the satellite elevation, so that we can visually inspect and
% assess directly how faithfully the main oscillation of each pass 
% is reproduced by the simulation (model accuracy)and identify systematic 
% differences due to multipath or model inaccuracies, we can explore how 
% parameters (eg antenna height, water level deviations) affect the interf. curve

figure('Name','Observed vs Simulated SNR (V/V) – Individual Arcs');

% --- Define the exact arcs you want to visualize ---
arcList = [2564, 1363, 591];    % <-- specify arc numbers directly
numArcs = numel(arcList);
nArcsAvailable = numel(snr);

% --- Precompute cumulative lengths for indexing simulation ---
arcLengths = cellfun(@numel, snr);
cumLen = [0; cumsum(arcLengths(:))];

tl = tiledlayout(numel(arcList), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
selected = 0;

for a = 1:numel(arcList)
    arcIdx = arcList(a);

    % --- Skip arcs with February data ---
    arc_time_vec = time{arcIdx};
    arc_dt = datetime(arc_time_vec, 'ConvertFrom', 'posixtime');
    if any(month(arc_dt) == 2)
        fprintf('Skipping arc %d (contains February data)\n', arcIdx);
        continue;
    end

    % --- Observed data ---
    obs_db = snr{arcIdx};                   
    obs_vv = 10.^(obs_db/20);               
    obs_vv_norm = obs_vv / max(obs_vv(~isnan(obs_vv)));
    obs_elev = elevation{arcIdx};

    % --- Corresponding simulated segment ---
    sim_start = cumLen(arcIdx) + 1;
    sim_end   = cumLen(arcIdx + 1);
    if sim_end > numel(result.phasor_interf)
        warning('Simulation shorter than expected for arc %d', arcIdx);
        continue;
    end

    sim_vv = abs(result.phasor_interf(sim_start:sim_end));
    sim_vv_norm = sim_vv / max(sim_vv(~isnan(sim_vv)));

    % --- Interpolation (if length mismatch) ---
    if numel(sim_vv_norm) ~= numel(obs_vv_norm)
        sim_idx = linspace(1, numel(sim_vv_norm), numel(sim_vv_norm));
        obs_idx = linspace(1, numel(sim_vv_norm), numel(obs_vv_norm));
        sim_vv_norm = interp1(sim_idx, sim_vv_norm, obs_idx, 'linear', 'extrap');
    end

    % --- Plot ---
    ax = nexttile;
    plot(ax, sin(deg2rad(obs_elev)), sim_vv_norm, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Simulated'); hold(ax, 'on');
    plot(ax, sin(deg2rad(obs_elev)), obs_vv_norm, 'b-', 'LineWidth', 1.2, 'DisplayName', 'Observed');
    title(ax, sprintf('Arc %d', arcIdx));
    grid(ax, 'on');
    if selected == 0
        legend(ax, 'Location', 'best');
    end
    selected = selected + 1;
end

% --- Shared axis labels ---
xlabel(tl, 'sin(Elevation)');
ylabel(tl, 'Normalized SNR (V/V)');
