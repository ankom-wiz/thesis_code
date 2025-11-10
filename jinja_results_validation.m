%% ======================================
% MPSim simulation for antenna heights and Water Level (WL) deviations
% Optimized to avoid memory issues
%% ======================================

%% Load GNSS merged arcs (observed data file)
load('merged_arcs.mat', 'snr', 'elevation', 'azimuth', 'time');

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

% Remove February (incomplete dataset)
isFeb = month(time_dt) == 2;
all_time_nofeb = all_time_sorted(~isFeb);
all_elev_nofeb = all_elev_sorted(~isFeb);
all_azim_nofeb = all_azim_sorted(~isFeb);
all_snr_nofeb  = all_snr_sorted(~isFeb);
time_dt_nofeb  = time_dt(~isFeb);

% Assign variables for simulation
time_ds  = all_time_nofeb;
elev_ds  = all_elev_nofeb;
azim_ds  = all_azim_nofeb;
snr_ds   = all_snr_nofeb;
time_dt_ds = time_dt_nofeb;

%% Load gauge data for water surface statistics
load('highResGaugeData.mat', 'highResGaugeData');

obsDates = dateshift(time_dt_ds, 'start', 'day');
gaugeDates = dateshift(highResGaugeData.Timestamp, 'start', 'day');

numObs = numel(obsDates);
wl_std_vec = nan(numObs,1);
uniqueDays = unique(obsDates);

for i = 1:numel(uniqueDays)
    modelIdx = find(obsDates == uniqueDays(i));
    gaugeIdx = find(gaugeDates == uniqueDays(i));
    if isempty(gaugeIdx), continue; end
    dev_day = highResGaugeData.Deviation(gaugeIdx);
    n_model = numel(modelIdx);
    n_dev   = numel(dev_day);
    chunk_sizes = floor(n_model/n_dev)*ones(n_dev,1);
    remainder = n_model - sum(chunk_sizes);
    chunk_sizes(1:remainder) = chunk_sizes(1:remainder)+1;
    idx_start = 1;
    for j = 1:n_dev
        idx_end = idx_start + chunk_sizes(j) - 1;
        wl_std_vec(modelIdx(idx_start:idx_end)) = dev_day(j);
        idx_start = idx_end + 1;
    end
end
wl_std_vec = abs(wl_std_vec);

%% Scenario definitions
antenna_heights = [2.6, 1.2, 10];  % meters
height_labels   = {'Baseline','Height – Low','Height – High'};
wldev_offsets = [0.3, 0.8];
wldev_labels  = {'WL deviation +0.3','WL deviation +0.8'};

%% Print header for table
fprintf('\n==================== SNR Validation Summary ====================\n');
fprintf('| %-15s | %-12s | %-12s | %-8s | %-4s | %-4s | %-5s |\n', ...
    'Scenario','Antenna H (m)','WL Dev (m)','Mean SNR','SD','R²','RMSE');
fprintf('-----------------------------------------------------------------\n');

%% -------------------
% 1) Antenna height simulations
for i = 1:numel(antenna_heights)
    h = antenna_heights(i);
    fprintf('\n=== Simulating antenna height: %.2f m ===\n', h);

    %% Setup MPSim
    sett = snr_settings();
    sett.ref.height_ant = h;
    sett.ant.slope = 'tipped';
    sett.ant.aspect = 176;
    sett.ant.model = 'isotropic horizontal';
    sett.ant.radome = 'NONE';
    sett.sfc.fnc_snr_setup_sfc_geometry = @snr_setup_sfc_geometry_horiz;
    sett.sfc.slope = 0;
    sett.sfc.aspect = 176;
    sett.sfc.dir_nrml = [0;0;1];
    sett.sfc.height_std = wl_std_vec;
    sett.sfc.vert_datum = 'top';
    sett.ref.pos_origin = [0,0,0];
    sett.sfc.pos_sfc0 = sett.ref.pos_origin;
    sett.sfc.fnc_snr_setup_sfc_material = @snr_setup_sfc_material_halfspaces;
    sett.sfc.material_top = 'air';
    sett.sfc.material_bottom = 'freshwater';
    sett.opt.gnss_name = 'GPS';
    sett.opt.freq_name = 'L1';
    sett.opt.rec.bandwidth_noise = 1;
    sett.opt.rec.temperature_noise = 300;
    sett.opt.rec.ant_density_noise_db = -199.3;
    sett.sat.num_obs = numel(elev_ds);
    sett.sat.elev = elev_ds(:);
    sett.sat.azim = azim_ds(:);
    sett.sat.epoch = time_ds(:);

    %% Run simulation
    setup = snr_setup(sett);
    [~, snr_db, ~, ~] = snr_fwd(setup);

    %% Compute statistics
    mean_snr = mean(snr_db, 'omitnan');
    sd_snr   = std(snr_db, 'omitnan');
    valid_idx = ~isnan(snr_ds) & ~isnan(snr_db);
    obs_valid = snr_ds(valid_idx);
    sim_valid = snr_db(valid_idx);
    RMSE     = sqrt(mean((sim_valid - obs_valid).^2));
    SS_res   = sum((obs_valid - sim_valid).^2);
    SS_tot   = sum((obs_valid - mean(obs_valid)).^2);
    R2       = 1 - SS_res / SS_tot;

    %% Print statistics immediately
    fprintf('| %-15s | %-12.2f | %-12s | %-8.2f | %-4.2f | %-4.2f | %-5.2f |\n', ...
        height_labels{i}, h, 'Gauge dev.', mean_snr, sd_snr, R2, RMSE);

    %% Clear memory for next scenario
    clear setup snr_db sim_valid obs_valid SS_res SS_tot;
    pack;  % optional: forces MATLAB to reclaim memory
end

%% -------------------
% 2) Surface WL deviation simulations (antenna fixed at 2.6 m)
for rIdx = 1:numel(wldev_offsets)
    offset = wldev_offsets(rIdx);
    label  = wldev_labels{rIdx};
    fprintf('\n=== Simulating WL deviation scenario: %s ===\n', label);

    %% Setup MPSim
    sett = snr_settings();
    sett.ref.height_ant = 2.6;
    sett.ant.slope = 'tipped';
    sett.ant.aspect = 176;
    sett.ant.model = 'isotropic horizontal';
    sett.ant.radome = 'NONE';
    sett.sfc.fnc_snr_setup_sfc_geometry = @snr_setup_sfc_geometry_horiz;
    sett.sfc.slope = 0;
    sett.sfc.aspect = 176;
    sett.sfc.dir_nrml = [0;0;1];
    sett.sfc.height_std = abs(wl_std_vec) + offset;
    sett.sfc.vert_datum = 'top';
    sett.ref.pos_origin = [0,0,0];
    sett.sfc.pos_sfc0 = sett.ref.pos_origin;
    sett.sfc.fnc_snr_setup_sfc_material = @snr_setup_sfc_material_halfspaces;
    sett.sfc.material_top = 'air';
    sett.sfc.material_bottom = 'freshwater';
    sett.opt.gnss_name = 'GPS';
    sett.opt.freq_name = 'L1';
    sett.opt.rec.bandwidth_noise = 1;
    sett.opt.rec.temperature_noise = 300;
    sett.opt.rec.ant_density_noise_db = -199.3;
    sett.sat.num_obs = numel(elev_ds);
    sett.sat.elev = elev_ds(:);
    sett.sat.azim = azim_ds(:);
    sett.sat.epoch = time_ds(:);

    %% Run simulation
    setup = snr_setup(sett);
    [~, snr_db, ~, ~] = snr_fwd(setup);

    %% Compute statistics
    mean_snr = mean(snr_db, 'omitnan');
    sd_snr   = std(snr_db, 'omitnan');
    valid_idx = ~isnan(snr_ds) & ~isnan(snr_db);
    obs_valid = snr_ds(valid_idx);
    sim_valid = snr_db(valid_idx);
    RMSE     = sqrt(mean((sim_valid - obs_valid).^2));
    SS_res   = sum((obs_valid - sim_valid).^2);
    SS_tot   = sum((obs_valid - mean(obs_valid)).^2);
    R2       = 1 - SS_res / SS_tot;

    %% Print statistics immediately
    fprintf('| %-15s | %-12.2f | %-12s | %-8.2f | %-4.2f | %-4.2f | %-5.2f |\n', ...
        label, 2.6, sprintf('+%.1f', offset), mean_snr, sd_snr, R2, RMSE);

    %% Clear memory for next scenario
    clear setup snr_db sim_valid obs_valid SS_res SS_tot;
end

fprintf('=================================================================\n');