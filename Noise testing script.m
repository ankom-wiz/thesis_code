%% MPSim  simulation (testing of noise parameters)
%% Load GNSS merged arcs (observed data file)
load('merged_arcs.mat', 'snr', 'elevation', 'azimuth', 'time');

%% Fix merged_mat parameter data forms
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

%% Downsampling option
ds = 1;
time_ds = all_time_nofeb(1:ds:end);
elev_ds = all_elev_nofeb(1:ds:end);
azim_ds = all_azim_nofeb(1:ds:end);
snr_ds  = all_snr_nofeb(1:ds:end);
time_dt_ds = time_dt_nofeb(1:ds:end);

%% Water surface statistics from high-resolution gauge data
load('highResGaugeData.mat', 'highResGaugeData');

obsDates = dateshift(time_dt_ds, 'start', 'day');
gaugeDates = dateshift(highResGaugeData.Timestamp, 'start', 'day');

numObs = numel(obsDates);
wl_std_vec = nan(numObs,1);
uniqueDays = unique(obsDates);

for i = 1:numel(uniqueDays)
    modelIdx = find(obsDates == uniqueDays(i));
    gaugeIdx = find(gaugeDates == uniqueDays(i));
    if isempty(gaugeIdx)
        continue;
    end
    
    dev_day = highResGaugeData.Deviation(gaugeIdx);
    n_model = numel(modelIdx);
    n_dev = numel(dev_day);
    chunk_sizes = floor(n_model / n_dev) * ones(n_dev,1);
    remainder = n_model - sum(chunk_sizes);
    chunk_sizes(1:remainder) = chunk_sizes(1:remainder) + 1;
    
    idx_start = 1;
    for j = 1:n_dev
        idx_end = idx_start + chunk_sizes(j) - 1;
        wl_std_vec(modelIdx(idx_start:idx_end)) = dev_day(j);
        idx_start = idx_end + 1;
    end
end

%% MPSim settings
sett = snr_settings();

%% Customize antenna and environment parameters
sett.ref.height_ant = 2.6;
sett.ant.slope = 'tipped';
sett.ant.aspect = 176;
sett.ant.model = 'isotropic horizontal';
sett.ant.radome = 'NONE';

%% Surface properties
sett.sfc.fnc_snr_setup_sfc_geometry = @snr_setup_sfc_geometry_horiz;
sett.sfc.slope = 0;
sett.sfc.aspect = 176;
sett.sfc.dir_nrml = [0;0;1];
sett.sfc.height_std = wl_std_vec;
sett.sfc.vert_datum = 'top';
sett.ref.pos_origin = [0, 0, 0];
sett.sfc.pos_sfc0 = sett.ref.pos_origin;
sett.sfc.fnc_snr_setup_sfc_material = @snr_setup_sfc_material_halfspaces;
sett.sfc.material_top = 'air';
sett.sfc.material_bottom = 'freshwater';

%% GNSS and signal configuration
sett.opt.gnss_name = 'GPS';
sett.opt.freq_name = 'L1';
sett.opt.rec.bandwidth_noise = 1;
sett.opt.rec.ant_density_noise_db = -199.3;

%% Observation directions
sett.sat.num_obs = numel(elev_ds);
sett.sat.elev = elev_ds(:);
sett.sat.azim = azim_ds(:);
sett.sat.epoch = time_ds(:);

%% Run 3 simulations with different environmental noise levels
noise_levels = [200, 300, 400];
labels = {'Cold receiver (200 K)', 'Nominal receiver (300 K)', 'Hot receiver (400 K)'};
colors = {'g', 'r', 'm'};  % for each simulation

figure('Name', 'Observed vs Simulated SNR (dB)');
days = dateshift(time_dt_ds,'start','day');
[unique_days, ~, idx_days] = unique(days);
daily_mean_obs = accumarray(idx_days, snr_ds, [], @mean);
plot(unique_days, daily_mean_obs, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Observed'); 
hold on;

for k = 1:length(noise_levels) 
    sett.opt.rec.temperature_noise = noise_levels(k);
    % Prepare and run simulation
    setup = snr_setup(sett);
    [~, snr_db, ~, ~] = snr_fwd(setup);
    
    % Compute daily mean
    daily_mean_sim = accumarray(idx_days, snr_db, [], @mean);
    
    % Plot
    plot(unique_days, daily_mean_sim, '-', 'Color', colors{k}, 'LineWidth', 1.5, ...
        'DisplayName', labels{k});
end

xlabel('Date');
ylabel('Daily Mean SNR (dB)');
title('Observed vs Simulated SNR under Different Noise Conditions');
legend('Location', 'best');
grid on;
