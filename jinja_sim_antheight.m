%% ======================================
% MPSim simulation for multiple antenna heights
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
% Simulation for Multiple Antenna Heights
antenna_heights = [1.2, 2.6, 10];  % meters
numHeights = numel(antenna_heights);

daily_mean_sim_all = nan(numel(unique(obsDates)), numHeights);
results_all = cell(numHeights,1);

for i = 1:numHeights
    h = antenna_heights(i);
    fprintf('\n=== Simulating antenna height: %.2f m ===\n', h);

    %% Step 1: Load default settings
    sett = snr_settings();

    %% Step 2: Customize antenna
    sett.ref.height_ant = h;
    sett.ant.slope = 'tipped';
    sett.ant.aspect = 176;
    sett.ant.model = 'isotropic horizontal';
    sett.ant.radome = 'NONE';

    %% Step 3: Surface properties
    sett.sfc.fnc_snr_setup_sfc_geometry = @snr_setup_sfc_geometry_horiz; 
    sett.sfc.slope = 0;
    sett.sfc.aspect = 176;
    sett.sfc.dir_nrml = [0;0;1];
    sett.sfc.height_std = wl_std_vec;
    sett.sfc.vert_datum = 'top';
    sett.ref.pos_origin = [0,0,0];
    sett.sfc.pos_sfc0 = sett.ref.pos_origin;

    % Reflecting surface
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
    results_all{i} = result;

    %% Step 7: Overall Statistics
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


    %% Step 8: Daily mean
    days = dateshift(time_dt_ds,'start','day');
    [unique_days, ~, idx_days] = unique(days);
    daily_mean_sim_all(:, i) = accumarray(idx_days, snr_db, [], @mean, NaN);
end

%% ======================================
% Plot Daily Mean SNR
figure('Name','Simulated SNR – Antenna Height Comparison');
shades = lines(numHeights);
hold on;
for i = 1:numHeights
    if antenna_heights(i) == 2.6
        plot(unique_days, daily_mean_sim_all(:, i), 'LineWidth', 1.8, 'Color', 'r', ...
        'DisplayName', sprintf('%.1f m', antenna_heights(i)));
    else
        plot(unique_days, daily_mean_sim_all(:, i), 'LineWidth', 1.2, 'Color', shades(i,:), ...
        'DisplayName', sprintf('%.1f m', antenna_heights(i)));
    end
end
xlabel('Date');
ylabel('Daily Mean SNR (dB)');
title('Simulated Daily Mean SNR for Multiple Antenna Heights');
legend('Location','best');
grid on;
hold off;

%% ======================================
% Individual Arc Analysis for Arc 591 – all antenna heights together
arcList = 591;  % Only Arc 591
arcLengths = cellfun(@numel, snr);
cumLen = [0; cumsum(arcLengths(:))];

figure('Name', 'Observed vs Simulated SNR – Arc 591, All Antenna Heights');
tl = tiledlayout(numHeights,1,'TileSpacing','compact','Padding','compact');  % one tile per antenna height

legendAdded = false;  % flag to add legend only once

for hIdx = 1:numHeights
    h = antenna_heights(hIdx);
    result_curr = results_all{hIdx};

    arcIdx = arcList;
    arc_time_vec = time{arcIdx};
    arc_dt = datetime(arc_time_vec, 'ConvertFrom', 'posixtime');
    
    if any(month(arc_dt) == 2)
        fprintf('Skipping arc %d (contains February data)\n', arcIdx);
        continue;
    end

    % Observed data
    obs_db = snr{arcIdx};
    obs_vv = 10.^(obs_db/20);
    obs_vv_norm = obs_vv / max(obs_vv(~isnan(obs_vv)));
    obs_elev = elevation{arcIdx};

    % Simulated segment
    sim_start = cumLen(arcIdx)+1;
    sim_end = cumLen(arcIdx+1);
    if sim_end > numel(result_curr.phasor_interf)
        warning('Simulation shorter than expected for arc %d (height %.1f m)', arcIdx, h);
        continue;
    end
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
    plot(ax, sin(deg2rad(obs_elev)), sim_vv_norm,'r-','LineWidth',1.5,'DisplayName','Simulated'); hold(ax,'on');
    plot(ax, sin(deg2rad(obs_elev)), obs_vv_norm,'b-','LineWidth',1.2,'DisplayName','Observed');
    title(ax,sprintf('Arc %d – Antenna Height %.1f m', arcIdx, h));
    grid(ax,'on');

    % Add legend only for the first tile
    if ~legendAdded
        legend(ax,'Location','best');
        legendAdded = true;
    end
end

xlabel(tl,'sin(Elevation)');
ylabel(tl,'Normalized SNR (V/V)');