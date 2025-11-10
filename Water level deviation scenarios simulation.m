%% ======================================
% MPSim simulation - Varying Water Level Deviation scenarios
%% Load GNSS merged arcs (observed data file)
load('merged_arcs.mat', 'snr', 'elevation', 'azimuth', 'time');

%% ======================================
% Fix merged_mat parameter data forms
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

% Remove February (incomplete dataset)
isFeb = month(time_dt) == 2;
all_time_nofeb = all_time_sorted(~isFeb);
all_elev_nofeb = all_elev_sorted(~isFeb);
all_azim_nofeb = all_azim_sorted(~isFeb);
all_snr_nofeb  = all_snr_sorted(~isFeb);
time_dt_nofeb  = time_dt(~isFeb);

%% ======================================
% Downsampling
ds = 1;   % downsampling factor
time_ds = all_time_nofeb(1:ds:end);
elev_ds = all_elev_nofeb(1:ds:end);
azim_ds = all_azim_nofeb(1:ds:end);
snr_ds  = all_snr_nofeb(1:ds:end);
time_dt_ds = time_dt_nofeb(1:ds:end);

%% ======================================
% Water surface statistics from high-resolution gauge data
load('highResGaugeData.mat', 'highResGaugeData');

obsDates = dateshift(time_dt_ds, 'start', 'day');
gaugeDates = dateshift(highResGaugeData.Timestamp, 'start', 'day');

numObs = numel(obsDates);
wl_std_vec = nan(numObs,1);
uniqueDays = unique(obsDates);

for iDay = 1:numel(uniqueDays)
    modelIdx = find(obsDates == uniqueDays(iDay));
    gaugeIdx = find(gaugeDates == uniqueDays(iDay));
    
    if isempty(gaugeIdx)
        continue;
    end
    
    dev_day = highResGaugeData.Deviation(gaugeIdx);
    n_model = numel(modelIdx);
    n_dev = numel(dev_day);
    
    % Divide model timestamps into equal chunks
    chunk_sizes = floor(n_model / n_dev) * ones(n_dev,1);
    remainder = n_model - sum(chunk_sizes);
    chunk_sizes(1:remainder) = chunk_sizes(1:remainder) + 1;
    
    % Assign deviations
    idx_start = 1;
    for j = 1:n_dev
        idx_end = idx_start + chunk_sizes(j) - 1;
        wl_std_vec(modelIdx(idx_start:idx_end)) = dev_day(j);
        idx_start = idx_end + 1;
    end
end

%% ======================================
% MPSim simulation for water level deviation scenarios
%% ======================================
% Define WL deviation scenarios
wldev_offsets = [0, 0.3, 0.8];   % meters
scenario_labels = {'+0.0m WL deviation', '+0.3m WL deviation', '+0.8m WL deviation'};
numScenarios = numel(wldev_offsets);

% Preallocate storage
daily_mean_sim_all = nan(numel(unique(obsDates)), numScenarios);
results_all = cell(numScenarios,1);

for rIdx = 1:numScenarios
    scenario_label = scenario_labels{rIdx};
    fprintf('\n=== Simulating WL deviation scenario: %s ===\n', scenario_label);

    %% Step 1: Load default settings
    sett = snr_settings();

    %% Step 2: Customize antenna
    sett.ref.height_ant = 2.6;
    sett.ant.slope = 'tipped';
    sett.ant.aspect = 176;
    sett.ant.model = 'isotropic horizontal';
    sett.ant.radome = 'NONE';

    %% Step 3: Surface properties
    sett.sfc.fnc_snr_setup_sfc_geometry = @snr_setup_sfc_geometry_horiz; 
    sett.sfc.slope = 0;
    sett.sfc.aspect = 176;
    sett.sfc.dir_nrml = [0;0;1];
    
    % Apply scenario-specific offset
    offset = wldev_offsets(rIdx);
    % Efficient offset application (no extra copies)
    sett.sfc.height_std = abs(wl_std_vec);                    % reference original vector
    sett.sfc.height_std = sett.sfc.height_std + offset;  % add scalar in-place
    
    sett.sfc.vert_datum = 'top';
    sett.ref.pos_origin = [0,0,0];
    sett.sfc.pos_sfc0 = sett.ref.pos_origin;

    sett.sfc.fnc_snr_setup_sfc_material = @snr_setup_sfc_material_halfspaces;
    sett.sfc.material_top = 'air';
    sett.sfc.material_bottom = 'freshwater';

    %% Step 4: GNSS configuration
    sett.opt.gnss_name = 'GPS';
    sett.opt.freq_name = 'L1';
    sett.opt.rec.bandwidth_noise = 1;
    sett.opt.rec.temperature_noise = 300;
    sett.opt.rec.ant_density_noise_db = -199.3;

    %% Step 5: Observation directions
    sett.sat.num_obs = numel(elev_ds);
    sett.sat.elev = elev_ds(:);
    sett.sat.azim = azim_ds(:);
    sett.sat.epoch = time_ds(:);

    %% Step 6: Setup and run simulation
    setup = snr_setup(sett);
    [result, snr_db, phasor_interf, phasor_composite] = snr_fwd(setup);
    results_all{rIdx} = result;

    %% Step 7: Overall statistics
    mean_snr = mean(snr_db, 'omitnan');
    median_snr = median(snr_db, 'omitnan');
    std_snr = std(snr_db, 'omitnan');
    min_snr = min(snr_db);
    max_snr = max(snr_db);
    p5_snr = prctile(snr_db,5);
    p25_snr = prctile(snr_db,25);
    p75_snr = prctile(snr_db,75);
    p95_snr = prctile(snr_db,95);

    fprintf('Mean: %.2f | Median: %.2f | SD: %.2f\n', mean_snr, median_snr, std_snr);
    fprintf('Min: %.2f | Max: %.2f\n', min_snr, max_snr);
    fprintf('5th: %.2f | 25th: %.2f | 75th: %.2f | 95th: %.2f\n', p5_snr, p25_snr, p75_snr, p95_snr);

    %% Step 8: Daily mean SNR
    days = dateshift(time_dt_ds,'start','day');
    [unique_days, ~, idx_days] = unique(days);
    daily_mean_sim_all(:, rIdx) = accumarray(idx_days, snr_db, [], @mean, NaN);
end

%% ======================================
% Plot Daily Mean SNR for WL Deviation Scenarios
figure('Name','Simulated SNR – WL Deviation Scenarios');
hold on;

colors = {'r','g','b'};  % Baseline, +0.3m, +0.8m WL Deviation
labels = {'+0.0m WL deviation', '+0.3m WL deviation', '+0.8m WL deviation'};

for rIdx = 1:numScenarios
    plot(unique_days, daily_mean_sim_all(:, rIdx), ...
        'LineWidth', 2.0, ...
        'Color', colors{rIdx}, ...
        'DisplayName', labels{rIdx});
end

xlabel('Date');
ylabel('Daily Mean SNR (dB)');
title('Simulated Daily Mean SNR – WL Deviation Scenarios');
legend('Location','best');
grid on;
hold off;

%% ======================================
% Individual Arc Analysis for Arc 591 – All WL Deviation Scenarios
arcIdx = 591;  % Only Arc 591
arcLengths = cellfun(@numel, snr);
cumLen = [0; cumsum(arcLengths(:))];

figure('Name', 'Observed vs Simulated SNR – Arc 591, All WL Deviation Scenarios');
tl = tiledlayout(numScenarios,1,'TileSpacing','compact','Padding','compact');
legendAdded = false;

arc_time_vec = time{arcIdx};
arc_dt = datetime(arc_time_vec, 'ConvertFrom', 'posixtime');

if any(month(arc_dt) == 2)
    fprintf('Skipping arc %d (contains February data)\n', arcIdx);
else
    % Observed data
    obs_db = snr{arcIdx};
    obs_vv = 10.^(obs_db/20);
    obs_vv_norm = obs_vv / max(obs_vv(~isnan(obs_vv)));
    obs_elev = elevation{arcIdx};

    for rIdx = 1:numScenarios
        result_curr = results_all{rIdx};

        sim_start = cumLen(arcIdx)+1;
        sim_end = cumLen(arcIdx+1);
        sim_vv = abs(result_curr.phasor_interf(sim_start:sim_end));
        sim_vv_norm = sim_vv / max(sim_vv(~isnan(sim_vv)));

        % Interpolate if mismatch
        if numel(sim_vv_norm) ~= numel(obs_vv_norm)
            sim_idx = linspace(1,numel(sim_vv_norm),numel(sim_vv_norm));
            obs_idx = linspace(1,numel(sim_vv_norm),numel(obs_vv_norm));
            sim_vv_norm = interp1(sim_idx, sim_vv_norm, obs_idx,'linear','extrap');
        end

        % Plot
        ax = nexttile;
        plot(ax, sin(deg2rad(obs_elev)), sim_vv_norm,'r-','LineWidth',1.2,'DisplayName','Simulated'); hold(ax,'on');
        plot(ax, sin(deg2rad(obs_elev)), obs_vv_norm,'b-','LineWidth',1.2,'DisplayName','Observed');
        title(ax, sprintf('Arc %d – %s', arcIdx, labels{rIdx}));
        grid(ax,'on');

        if ~legendAdded
            legend(ax,'Location','best');
            legendAdded = true;
        end
    end
end

xlabel(tl,'sin(Elevation)');
ylabel(tl,'Normalized SNR (V/V)');
