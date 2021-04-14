function spikenet50_FICspike_trig()
% Runs an integerate and fire neural network model.
% Displays correlograms and spike triggered averages for selected units.
%
% This version has been setup for spike triggered conditioning with 
% Fixed Intracolumn Connections (FIC).  FIC assumes the role of the
% correlated bias inputs to provide correlated activity within
% each column.  This is a modification of spikenet50.m with:
%   p.conditioning_type = 1;
%   p.bias_corrrelated_fraction = 0.0
% Full connectivity for intracolumn connections fixed between 200-300 uV.
%   p.initmin = 200;
%   p.initmax = 300;
%   p.epsp2excit_incol = [p.initmin p.initmax 1 0];     % Epsp for connections to excitatory units within a column
%   p.epsp2inhib_incol = [p.initmin p.initmax 1 0];     % Epsp for connections to inhibitory units within a column
%   p.ipsp2excit_incol = [-p.initmax -p.initmin 1 0];   % Ipsp for connections to excitatory units within a column
%   p.ipsp2inhib_incol = [-p.initmax -p.initmin 1 0];   % Ipsp for connections to inhibitory units within a column
% Stimulation parameters:
%   p.stim_delay_ms = 10;    % Millisecond delay between spike on Ae1 and stimulus on Column B.
%   p.stim_pulse_train = 1;  % Can be 2 or 3 for stimulus trains of 2 or 3 pulses per train.
%   p.stim_uV = 2000;        % Size of conditioning stimulus (0 can be used as a sham stimulus)
%   p.stim_refractory_ms = 10;  % Refractory period on delivered stimulation.
%
%-----------------------------
%
% The model is comprised of integrate-and-fire units divided into groups 
% referred to as Columns A, B, and C.  Units receive excitatory background
% activity proportioned between uncorrelated and within-column correlated 
% events.  Connections from source to target units are assigned sparsely
% with inhibitory units limited to within-column targets.  Each unit
% maintains a potential which represents the sum of the synaptic inputs to 
% that unit.  When a source unit’s potential exceeds a threshold, that
% unit “fires” resetting its potential to zero and initiating a spike, 
% which arrives after a conduction delay, to all its target units.  Each
% spike arriving at a unit will add a weighted function to that unit’s 
% potential modeling the post-synaptic potentials (PSPs) of physiological
% neurons.  When spike-timing dependent plasticity is active, connection
% weights are modified based on the differences between target and source
% unit firing times. Networks can be run under various conditioning methods
% to simulate experimental stimulation protocols. 
%
%-----------------------------

close all;     % Can be used to close all currently open figures.
[p.scriptpath, p.scriptname] = fileparts(mfilename('fullpath'));
cd(p.scriptpath); % Run inside of this scripts directory.

% Output is saved to a folder named [p.folder '\' p.subfolder '\' p.prefix '_' datestring '_' indexstring]
p.folder = 'c:\data\spikenet'; % Data directory prefix
p.subfolder = 'spikenet50';       % Folder name to hold ouput folders.
p.prefix = 'sn';               % Output folder name prefix.
p.mexname = 'spikenet50mex';   % Name of the mex routine that runs the network.

% Simulation parameters.  Note that not all parameters affect all conditioning types.
p.conditioning_type = 1;       % 0=No conditioning, 1=spike triggered, 2=paired pulse, 3 = cycle triggered on LFPB, 4 = Tetanic stim on Col A, 5 = Tetanic stim on Col B, 6 = Tetatic stim on Col C, 7 = EMG triggered Stimulation, 8 = Gamma triggered stimulation, 9 = Activity triggered stimulation
p.stim_delay_ms = 10;          % Stimulaton delay from rec_unit threshold crossing to start of stimulation PSP in designated stim target units.
p.stim_pulse_train = 1;        % Pulses in the 300 Hz stim pulse train for conditioning. Limited to 3 for Paired Pulse stimulation.
p.stim_pulse_freq = 300;       % Interstimulus interval for stimlus trains (this will be converted to an integer number of timesteps)
p.stim_phase = 0;              % Cycle triggered stim phase, 0 = rising, 90 = peak, 180 = falling, 270 = trough. These phases are fairly exact, others are interpolated from these points. Phases between 46..89 and 316..359 have the most interpolation variability.
p.stim_refractory_ms = 10;     % Refractory period on spike triggered and LFP cycle triggered stimulation
p.lfp_detect_level = 30000;    % Amplitude of LFP used for cycle triggered conditioning (use 30000 for cycle triggered +/-20%, -20000 or so for gamma triggered depending on %correlated bias drive, 1000 for EMG triggered
p.gamma_band = [50 80];        % Bandwidth for conditioning_type in Hz (for example [50 80] for 50Hz to 80Hz bandpass filter.
p.emg_band = [100 2500];       % Column A motor output filter (Hz) for EMG triggered stimulation.
p.conditioning_secs = 10;      % Seconds in each trial used for conditioning.  Limits conditioning to beginning of each 10 second trial.
p.bias_modulation_rate_A = 0;  % Cycles per second modulation rate on Column A. 0 = no modulation. 19, 20, or 21 for sine wave
p.bias_modulation_rate_B = 0;  % Cycles per second modulation rate on Column B. 0 = no modulation. 20 for 20Hz sine wave
p.beta_band = [15 25];         % LFP_B filter band (Hz) for cycle triggered stimulation.
p.bias_modulation_rate_C = 0;  % Cycles per second modulation rate on Column C. 0 = no modulation. 19, 20 or 21 for sine wave
p.bias_modulation_amount = 0.0;  % Size of modulation. 0.2 = +/- 20% of normal rate.  Normally 0. Cycle triggered stim uses 0.2, or set to 0 to just have test pulses for comparison purposes.  Activity dependent stim usues 0.4, and tetanic will use anything as a control for other conditioning types
p.bias_modulation_step_secs = 2; % Number of seconds between modulation episodes.  Normally 2, use 5 to test for training decay.
p.tetanic_freq = 10;           % Stims per sec for tetanic and activity dependent conditioning, Max stims/sec for exponentially distributed tetanic conditioning (actual rate is always lower because of refractory period)
p.test_pulse_train = 1;        % Number of pulses in the stimulus test train
p.test_pulse_freq = 300;       % Test pulse rate in Hz (300Hz gives 3.3 ms interval, 500Hz gives 2 ms interval)

p.test_A1_lesion = 0;          % 1 if testing should be done with weights from unit A1 set to 0.
p.test_AB_lesion = 0;          % 1 if testing should be done with weights from Col A to Col B set to 0.
p.test_AC_lesion = 0;          % 1 if testing should be done with weights from Col A to Col C set to 0.
p.test_No_ColA_ColB_connections = 0; % 0 Normal connectivity, 1 = disallow all connections between Col A and ColB.
p.do_raster_figure = 0;        % 1 to calculate LFP power spectrum and show along with a dot raster.  This can be slow.

% File for replaying stimulus events.  Empty string for no replay.
% This can be used for spike triggered and cycle triggered catch trials.
%p.stim_replay_file = 'C:\data\spikenet\spikenet37\sn37_20180209_02\stim_output_clockticks.mat';
p.stim_replay_file = '';

% Number of trials to run under various conditions.
p.ntrials_off = 50;     % Number of initial trials with conditioning off (STDP training on).
p.ntrials_on = 50;      % Number of following trials with conditioning on (STDP training on).
p.ntrials_test = 50;    % Number of trials for summary figures (both conditioning and STDP training off)
p.ntrials_repeat = 1;   % Number of repetitions through the off/test/on/test cycle.

% Network configuration constants.  At this point, these cannot be changed
% without making changes to most of the analysis and figures.
p.n_excit = 40;  % Number of units in each excitatory group
p.n_inhib = 40;  % Number of units in each inhibitory group.
p.n_out = 40;    % Number of final output units for each Column
p.n_cols = 3;    % Number of simulated columns.

% Threshold for cortical and motor output units.
p.threshold = 5000;        % uV threshold is currently global for all cortical units
p.output_threshold_range = [5000 6000]; % Output units have graded thresholds to simulate connections to muscles.
p.output_field_range = [500 1500];      % Output unit graded strengths to the output EMG.

% Axonal and dendritic delays affect the delivery time of post-synaptic
% poentials and the Spike Time Dependent Plasticity calculations.
% delay = axonal + dendritic delay.
%
% For any spike occurring on unit i at time t
%    For all weights Wij     (from unit j to unit i)
%      At time t:
%        1) The weakening potential Ti is updated by the weakening factor (c*r).
%        2) The strengthinging STDP rule is applied to Wij using the strengthing potential Sj
%    For all weights Wki     (from unit i to unit k)
%      At time t + delay:
%        1) Weight Wki is applied to Vk to initiate a PSP on unit k.
%        2) The strengthing potential Si is updated by the strengthing factor (r).
%        3) The weakinging STDP rule is applied to all Wki using the weakening potential Tk
%
% Total conduction delay is limited to 10 milliseconds. This could be increased to 20 ms
% by doubling the size of the spike_delay_queue in the spikenet mex function.

p.axonal_delay_ms = 2;     % Axonal mS delay is currently global for all non-output connections.
p.dendritic_delay_ms = 1;  % Dendritic mS delay is curretly global for all non-output connections
p.output_delay_ms = 10;    % Connection delay for connections to output units 

% PSP shape is exp(-t/tau_slow) - exp(-t/tau_fast). This will be converted
% into a pair of leaky integrators with Euler's method for fast computation
% using a timestep of 1/p.fs seconds

p.fs = 10000;  % Time steps per second. Keep as a multiple of 10000 with a maximum of 100000
p.PSP_slow_time_constant_ms = 3.2; % Default synaptic potential shape
p.PSP_fast_time_constant_ms = 0.8; % Time constants are in milliseconds

% Bias activity

p.bias_value = 350;                % uV strenth of bias potientials
p.bias_rate = 1800;                % Bias spikes per second to each column unit. Consider smaller values if biases are modulated.
p.bias_corrrelated_fraction = 0.0; % Fraction of correlated bias spikes
p.correlated_bias_std_ms = 3;      % <= 5ms. Standard deviation of the normally distributed (correlated) bias spikes (in milliseconds)
p.out_bias_rate = 2000;            % Number of uncorrelated bias spikes per second for output units.

p.correlated_bias_rate = p.bias_corrrelated_fraction * p.bias_rate; % Correlated bias spikes per second delivered to each unit in a column. These will cause common input peaks in cross-correlations of same column units.
p.uncorrelated_bias_rate = p.bias_rate - p.correlated_bias_rate;    % Number of bias spikes per second for normal firing rate.

% Weight limits.  Weight dependence slows down weight changes as a weight
% approaches zero or maximum (or -maximum for inhibitory weights).  This is
% an experimental way to apply some homeostasis to the weights to help keep
% them more central. Standard simulations use no weight dependence.
% If weight_dependence is non-zero, strengthening changes are multipled by
% (1 - (current_weight / weight_limit)) ^ weight_dependence, and weakening
% changes are multipled by (current_weight / weight_limit) ^ weight_dependence

p.max_strength = 500;     % Maximum allowed connection strength (uV)
p.weight_dependence = 0.0; % 0.001 for almost no weight dependence. Maximum of 1 for linear weight dependence.

% Initial PSP ranges in uV.  Arrays are [uVmin uVmax pChanceForConnection enableTraining]
p.initmin = 200;
p.initmax = 300;
p.epsp2excit_incol = [p.initmin p.initmax 1 0];     % Epsp for connections to excitatory units within a column
p.epsp2inhib_incol = [p.initmin p.initmax 1 0];     % Epsp for connections to inhibitory units within a column
p.epsp2excit_outcol = [p.initmin p.initmax 1/6 1];    % Epsp for connections to excitatory units in adjacent columns
p.epsp2inhib_outcol = [p.initmin p.initmax 1/6 1];    % Epsp for connections to inhibitory units in adjacent columns
p.ipsp2excit_incol = [-p.initmax -p.initmin 1 0];   % Ipsp for connections to excitatory units within a column
p.ipsp2inhib_incol = [-p.initmax -p.initmin 1 0];   % Ipsp for connections to inhibitory units within a column
p.epsp2output_incol = [350 350 1/3 0];    % Epsp for connections to output units.  These are not usually trained. Negative connection chance means first p * N units rather than random chance.

p.training_rate = 100;    % Train rate factor for both strengthing (pos) and weakening (neg) sides of the SDTP rule.
p.train_weakening = 0.55; % Relative amplitude for the weakening side of the SDTP curve.
p.STDP_strengthening_slow_time_constant_ms = 15.4; % STDP strengthening potential shape
p.STDP_strengthening_fast_time_constant_ms = 2;
p.STDP_weakening_slow_time_constant_ms = 33.3;     % STDP weakening potential shape
p.STDP_weakening_fast_time_constant_ms = 2;

% Conditioning related parameters

p.stim_uV = 2000;        % Size of conditioning stimulus (0 can be used as a sham stimulus)
if p.conditioning_type == 2
    p.pair_uV = 2000;    % Size of the paired pulse simulation (the second pulse, this is only for paired pulse stimulation)
else
    p.pair_uV = p.stim_uV;
end
p.test_uV = 3000;     % Size of test stimulus

stim_excit_range = (1:p.n_excit);   % Stimulate all excitatory units in target column
%stim_excit_range = []; % Use this line to disable stimulus on excitatory units.
stim_inhib_range = (1:p.n_inhib);   % Stimulate all inhibitory units in target column
%stim_inhib_range = []; % Use this line to disable stimulus on inhibitory units.

if p.conditioning_type == 1  % Spike Triggered Stimulation
    p.trigger_unit = 1;      % Source unit index for spike triggered stimulation.  0 will disable spike triggered stimulation.
    %p.trigger_unit = p.n_excit + p.n_inhib + 1;  % This line tests triggering from the first AOutput unit
else
    p.trigger_unit = 0;
end

% Define stimulation sources for non-delayed stimulation when using
% Tetanic, Cycle-Triggered, Activity-Triggered, and Paired pulse (first
% pulse) conditioning
if (p.conditioning_type == 2) || (p.conditioning_type == 3) || (p.conditioning_type == 4) || (p.conditioning_type == 9)
    p.stim_source_def = {'Ae'; stim_excit_range; 'Ai'; stim_inhib_range};  % Tetanic A, Paired Pulse, Activity triggered, or Cycle triggered stim 
elseif (p.conditioning_type == 5) 
    p.stim_source_def = {'Be'; stim_excit_range; 'Bi'; stim_inhib_range};  % Tetanic B stim 
elseif (p.conditioning_type == 6)
    p.stim_source_def = {'Ce'; stim_excit_range; 'Ci'; stim_inhib_range};  % Tetanic C stim 
else
    p.stim_source_def = {}; % No sources to stimulate
end

% Define stimulation targets for delay stimulation when using Spike- Cycle-
% EMG- Gamma- Triggered, and Paired Pulse (second pulse) conditioning
if (p.conditioning_type == 1) || (p.conditioning_type == 2) || (p.conditioning_type == 7) || (p.conditioning_type == 8)
    p.stim_target_def = {'Be'; stim_excit_range; 'Bi'; stim_inhib_range};
elseif (p.conditioning_type == 3)
    % Cycle triggered stimulation when using LFPB as trigger and Col A as stimulation target.   
    p.stim_target_def = {'Ae'; stim_excit_range; 'Ai'; stim_inhib_range};
else
    p.stim_target_def = {}; % No targets to stimulate
end

p.random_seed = 1; % Restart random number generator so similar networks topologies will have same layout.

%----------------------------- 
% In general, users will not change anything below this line unless adding
% or changing features.

p.start_time = clock();    % Mark start time.
rng(p.random_seed);        % Seed random number generator
p.version_id = mfilename;  % Use file name as a program identifier.
p.fs = min(100000, max(10000, round(p.fs / 10000) * 10000)); % Make sure p.fs is multiple of 10000 between 10000 and 100000.
p.msfs = round(p.fs / 1000); % Time steps per millisecond.
p.axonal_delay = floor(p.axonal_delay_ms * p.msfs);  % axonal psp delay in timesteps
p.dendritic_delay = floor(p.dendritic_delay_ms * p.msfs);  % dendridic psp delay in timesteps
p.conduction_delay = p.axonal_delay + p.dendritic_delay;   % total conduction time in timesteps
p.output_delay = floor(p.output_delay_ms * p.msfs);        % Connection delay for any connection to an output unit.
p.stim_delay = floor(p.stim_delay_ms * p.msfs);            % stimulation delay in timesteps
p.stim_refractory = floor(p.stim_refractory_ms * p.msfs);  % stimulation refractory
[p.gamma_filter_b, p.gamma_filter_a] = butter(1, p.gamma_band / (p.fs / 2)); % LFP_A filter for gamma triggered stimulation
[p.beta_filter_b, p.beta_filter_a] = butter(1, p.beta_band / (p.fs / 2)); % LFP_B filter for cycle triggered stimulation
[p.emg_filter_b, p.emg_filter_a] = butter(1, p.emg_band / (p.fs / 2)); % Motor output filter for emg triggered stimulation

% Convert time constants into decay factors for PSP and STDP shapes
% (this is Eulers method used for estimating an exponential decay function)
p.timestep2ms = 1000 / p.fs; % Milliseconds per timestep
p.psp_slow_decay = 1 - p.timestep2ms / p.PSP_slow_time_constant_ms; 
p.psp_fast_decay = 1 - p.timestep2ms / p.PSP_fast_time_constant_ms;   
p.train_pos_slow_decay = 1 - p.timestep2ms / p.STDP_strengthening_slow_time_constant_ms;   % Shape of strengthening STDP rule.
p.train_pos_fast_decay = 1 - p.timestep2ms / p.STDP_strengthening_fast_time_constant_ms;
p.train_neg_slow_decay = 1 - p.timestep2ms / p.STDP_weakening_slow_time_constant_ms;       % Shape of weakening STDP rule.
p.train_neg_fast_decay = 1 - p.timestep2ms / p.STDP_weakening_fast_time_constant_ms;

p.train_pos_factor = p.training_rate; % STDP Strengthening factor (r) for spike pairs where (tPost - tPre) >= 0.
p.train_neg_factor = p.train_weakening * p.train_pos_factor; % STDP weakening factor (cr) for spike pairs where (tPost - tPre < 0).

% Generate a normal distribution for the correlated spikes.  We clip this
% distribution at +/- 4 standard deviations and then center it at 
% 4 standard deviations.  Values in the table are then in timesteps.
% When a correlated input bias occurs for a column, all units in that
% column will receive a bias spike at an independent time taken from 
% this table.
p.correlated_bias_max_timesteps = 8 * p.correlated_bias_std_ms * p.msfs;
p.normal_pdf_table = round(random('norm', 0.5 * p.correlated_bias_max_timesteps, p.correlated_bias_std_ms * p.msfs, [100000, 1])); % used for correlated bias spikes
p.normal_pdf_table((p.normal_pdf_table < 0) | (p.normal_pdf_table > p.correlated_bias_max_timesteps)) = [];  % remove elements more than 4 standard deviations from mean.
p.normal_pdf_table = uint16(p.normal_pdf_table(1:2^16));
coffset = round(mean(p.normal_pdf_table));  % clocktick offset for any modulation of correlated bias inputs.  This is just the mean all possible time adjustments.

% Define the number of biases, units, and weights.
% Biases are used to give units background activity.  This activity can be
% modulated by the bias_chance array assigned to each bias.
% Steps is the number of time steps to run the network on each iteration.
% This value is usually in the 1 to 10 second range, but may be set to
% a time based on the length of an experimental trial in a simulated task.

p.units = p.n_cols * (p.n_excit + p.n_inhib + p.n_out);  % Number of actual units, but will be adjusted up if necessary.
p.weights = 50000;        % Initial number of weights, but will be adjusted up or down if necessary.
p.time_steps = floor(10 * p.fs); % 10 Seconds of simulation time per iteration.
p.train_info = zeros(1,4);% Place for keeping track of some training summary results

% Allocate space for biases.  Biases return a randomized value on every
% access using a 32-bit random number generator R.
% if R(access) < bias_chance(step), return bias_strength, else return 0.
% Bias 1 is for uncorrelated biases. Biases 2,3,4 are for correlated biases
% assigned to groups A,B,C.  Bias 5 is for output units.

p.bias_strength = zeros(5, 1, 'double'); % 
p.bias_chance = zeros(5, p.time_steps, 'uint32'); % probability(step) = chance(step)/(2^32 - 1)

% Pre-allocate space for units

p.unit_names = repmat({}, p.units, 1);            % Text name of unit
p.unit_group_index = zeros(p.units, 1, 'uint16'); % Unit Index number within its group
p.unit_threshold = zeros(p.units, 1, 'double');   % Each unit has its own threshold
p.unit_bias_offset = zeros(p.units, 1, 'uint16'); % zero based bias id
p.unit_pre_count = zeros(p.units, 1, 'uint16');   % Number of presynaptic units
p.unit_post_count = zeros(p.units, 1, 'uint16');  % Number of postsynaptic units
p.unit_pre_offset = zeros(p.units, 1, 'uint32');  % Zero based offset into weight_pre_sort
p.unit_post_offset = zeros(p.units, 1, 'uint32'); % Zero based offset into wieght_post_sort
p.unit_lfp_offset = zeros(p.units, 1, 'uint16');  % Zero based offset into LFP array
p.unit_column_id = zeros(p.units, 1, 'uint16');   % 0=Col A unit, 1 = B, 2 = C, output units have 16384 added
p.unit_stim_source = zeros(p.units, 1, 'uint16'); % Flag for each unit. 1 if used as a stim source.  Filled in by InitDone()
p.unit_stim_target = zeros(p.units, 1, 'uint16'); % Flag for each unit. 1 if used as a stim target.  Filled in by InitDone()
p.unit_output_field_strength = zeros(p.units, 1, 'double'); % Strength for EMG output.

% Pre-allocate space for weights

p.weight_pre_unit = zeros(p.weights, 1, 'uint16');      % Index of presynaptic unit (1..p.units)
p.weight_post_unit = zeros(p.weights, 1, 'uint16');     % Index of postsynpatic unit
p.weight_strength = zeros(p.weights, 1, 'double');      % Current connection strength
p.weight_training_rule = zeros(p.weights, 1, 'uint16'); % Learning rule type (0 = none);
p.weight_test_lesion = zeros(p.weights, 1, 'uint16');   % Flag for each weight.  1 if weight should bet set to 0 while testing.

% p.stim_source and target_times are used for paired pulse stimulations and tetanic stimulation. 
% p.stim_test_times are used for testing evoked potentials with a column
% wide stimulus during test sections.

p.stim_source_times = zeros(p.time_steps, 1, 'uint16');  % Default to no paired pulse or cycle triggered stims
p.stim_target_times = zeros(p.time_steps, 1, 'uint16');
p.stim_test_times = zeros(p.time_steps, 1, 'uint16');
p.stim_output_clockticks = zeros(0, 1, 'int32');         % Collect clocktick of each output stimulation.
p.unit_spike_counts = zeros(1, p.units);                 % collect number of  output spikes for each unit
p.stim_pulse_isi = floor(p.fs / p.stim_pulse_freq);      % stimulus pulse train interval.

test_pulse_isi = floor(p.fs / p.test_pulse_freq); % Test pulse ISI in timesteps
for istim = 0:p.test_pulse_train-1
    p.stim_test_times(floor(8.0 * p.fs + 1 + test_pulse_isi * istim)) = 1; % Col A test stimulus times start at 8 seconds
    p.stim_test_times(floor(8.7 * p.fs + 1 + test_pulse_isi * istim)) = 2; % Col B test stimulus times start at 8.7 seconds
    p.stim_test_times(floor(9.4 * p.fs + 1 + test_pulse_isi * istim)) = 3; % Col C test stimulus times start at 9.4 seconds
end

if (p.conditioning_type == 2)  %  paired pulse stim times.
    for episode = 1:p.conditioning_secs  % These stims need to fall outside of any modulation episodes
        offset = floor((episode - 1) * p.fs + 1);        
        p.stim_source_times(floor(offset + 0.1 * p.fs)) = 1;  % Stim Column A, Paired pulse is two sets of pairs per second
        p.stim_target_times(floor(offset + 0.1 * p.fs + p.stim_delay)) = 1;  % Paired pulse Delayed stim on Column B after each stim on Column A
        p.stim_source_times(floor(offset + 0.3 * p.fs)) = 1;
        p.stim_target_times(floor(offset + 0.3 * p.fs + p.stim_delay)) = 1;
        if p.stim_pulse_train >= 2
            p.stim_source_times(floor(offset + 0.1 * p.fs + p.stim_pulse_isi)) = 1;  % Stim Column A, Paired pulse is two sets of pairs per second
            p.stim_target_times(floor(offset + 0.1 * p.fs + p.stim_pulse_isi + p.stim_delay)) = 1;  % Paired pulse Delayed stim on Column B after each stim on Column A
            p.stim_source_times(floor(offset + 0.3 * p.fs + p.stim_pulse_isi)) = 1;
            p.stim_target_times(floor(offset + 0.3 * p.fs + p.stim_pulse_isi + p.stim_delay)) = 1;
        end
        if p.stim_pulse_train >= 3
            p.stim_source_times(floor(offset + 0.1 * p.fs + 2*p.stim_pulse_isi)) = 1;  % Stim Column A, Paired pulse is two sets of pairs per second
            p.stim_target_times(floor(offset + 0.1 * p.fs + 2*p.stim_pulse_isi + p.stim_delay)) = 1;  % Paired pulse Delayed stim on Column B after each stim on Column A
            p.stim_source_times(floor(offset + 0.3 * p.fs + 2*p.stim_pulse_isi)) = 1;
            p.stim_target_times(floor(offset + 0.3 * p.fs + 2*p.stim_pulse_isi + p.stim_delay)) = 1;
        end
    end
end

% Possibly load a file of recorded stimulation times to play back.
p.stim_replay_clockticks = -1;  % Last value in array must be < 0
if ~isempty(p.stim_replay_file)
    try
        S = load(p.stim_replay_file, 'clockticks');
        p.stim_replay_clockticks = S.clockticks - 1; % Convert to zero based indexing for the mex file.
        p.stim_replay_clockticks(end+1) = -1; % terminate the list
    catch
        p.stim_replay_clockticks = -1;
        disp(['Could not replay conditioning stimulations from file: ' p.stim_replay_file]);
    end
end

% Initialization

Init_Network();

%%%%
% Fill in bias information
% Columns can have their own bias activity profiles.
% Initialize with same correlated and uncorrelated rates for each column.

uncorrelated_rate1 = repmat(p.uncorrelated_bias_rate, p.time_steps, 1);
uncorrelated_rate2 = uncorrelated_rate1;
uncorrelated_rate3 = uncorrelated_rate1;
correlated_rate1 = repmat(p.correlated_bias_rate, p.time_steps, 1);
correlated_rate2 = correlated_rate1;
correlated_rate3 = correlated_rate1;
out_rate = repmat(p.out_bias_rate, p.time_steps, 1); % Output base bias rate

% Base_rate can be modulated here if cyclical activity is needed.
% This iscycle triggered conditioning and activity dependent stimulaton, but other conditioning
% types can have the modulations as well. Cyclical episodes can start at
% 0.5, 3.5, 7.5 seconds and have 7 oscillations each.

if (p.conditioning_type == 9) || (p.conditioning_type == 4) || (p.conditioning_type == 5) || (p.conditioning_type == 6)
    % Activity triggered stimulation pairs increased input bias on
    % Column B with 1 second of stimulation at p.tetanic_freq on Column A.
    % Tetanic conditioning is allowed to have the modulated activity
    % increase on Column B for comparisons.
    ramp = [(0:0.2*p.fs-1)/(0.2 * p.fs), ones(1, p.fs), (0.2*p.fs-1:-1:0)/(0.2 * p.fs)] * p.bias_modulation_amount + 1;
    stims = [];
    if p.tetanic_freq == 1
        stims = 0;
    elseif p.tetanic_freq > 0   
        stims = floor((0:(p.tetanic_freq-1)) * p.fs / p.tetanic_freq);
        % The following lines are for testing short 0.2 seconds of stim
        %stims = 0:floor(p.fs / p.tetanic_freq):0.2*p.fs;
        %stims = stims + floor(0.5 * (0.2*p.fs - stims(end))); % Center the train
    end
    for iramp = 0:2
        % Place 3 of these ramps in our 10 second trial time.
        ioffset = floor(iramp * 2 * p.fs + 1);  % Ramp runs every 2 seconds
        range = ioffset:(ioffset+length(ramp)-1);
        correlated_rate2(range) = correlated_rate2(range) .* ramp';
        if ~isempty(stims) && (p.conditioning_type == 9) % Only setup stimulations for Activity Triggered conditioning
        	p.stim_source_times(ioffset + stims + p.stim_delay) = 1; % p.stim_delay needs to account for the 200ms ramp
        end
    end 
end

if p.bias_modulation_rate_A > 0
    sine_mod = 0.5 * (1 + sin((1:p.time_steps) * 2 * pi * p.bias_modulation_rate_A / p.fs));
    sine_mod = (1 - p.bias_modulation_amount) + 2 * p.bias_modulation_amount * sine_mod';  % Modulate between 80% to 120% of base uncorrelated rate.
    for episode = 1:p.bias_modulation_step_secs:p.conditioning_secs
        for cycle_start = floor(0.5 * p.fs)  % 1 set of 7 oscillations starting .5 seconds into each episode
            offset = floor((episode - 1) * p.fs + 1);
            index1 = offset + cycle_start;
            index2 = index1 + floor(7 * p.fs / p.bias_modulation_rate_A);
            uncorrelated_rate1(index1:index2) = uncorrelated_rate1(index1:index2) .* sine_mod(1:index2-index1+1);
            correlated_rate1(index1-coffset:index2-coffset) = correlated_rate1(index1-coffset:index2-coffset) .* sine_mod(1:index2-index1+1);
        end
    end
end

p.postlfp_test_sample = floor(0.25 * p.fs + 1);  % Conditioning type 3 cycle triggered pre-oscillation test pulse time.
p.prelfp_test_sample = floor(p.conditioning_secs * p.fs + .1 * p.fs + 1);  % Conditioning type 3 cycle triggered post-oscillation test pulse time.
p.test_index = floor([.3 1.1 2.3] * p.fs); % Test pulses for cycle triggered conditioning, [pre-cycles test, post-cycles test, next pre-cycle test]
if p.bias_modulation_rate_B > 0
    sine_mod = 0.5 * (1 + sin((1:p.time_steps) * 2 * pi * p.bias_modulation_rate_B / p.fs));
    sine_mod = (1 - p.bias_modulation_amount) + 2 * p.bias_modulation_amount * sine_mod';  % Modulate between 80% to 120% of base uncorrelated rate.
    p.stim_phase_offset = floor((p.fs / p.bias_modulation_rate_B) * p.stim_phase / 360); % 20Hz with 0.1ms time step is 500 steps in 360 degrees.
    for episode = 1:p.bias_modulation_step_secs:p.conditioning_secs-p.bias_modulation_step_secs
        for cycle_start = floor(0.5 * p.fs)
            offset = floor((episode - 1) * p.fs + 1);
            index1 = offset + cycle_start;
            index2 = index1 + floor(7 * p.fs / p.bias_modulation_rate_B);
            uncorrelated_rate2(index1:index2) = uncorrelated_rate2(index1:index2) .* sine_mod(1:index2-index1+1);
            correlated_rate2(index1-coffset:index2-coffset) = correlated_rate2(index1-coffset:index2-coffset) .* sine_mod(1:index2-index1+1);
            if (p.conditioning_type == 3) && (p.bias_modulation_step_secs < 5)
                p.stim_target_times(floor(p.test_index + offset)) = 1;  % Pre and post cycle test pulses go into target stim_target_time so they will not be marked as conditioning triggers.
            end
        end
    end
end

if p.bias_modulation_rate_C > 0
    sine_mod = 0.5 * (1 + sin((1:p.time_steps) * 2 * pi * p.bias_modulation_rate_C / p.fs));
    sine_mod = (1 - p.bias_modulation_amount) + 2 * p.bias_modulation_amount * sine_mod';  % Modulate between 80% to 120% of base uncorrelated rate.
    for episode = 1:p.bias_modulation_step_secs:p.conditioning_secs
        for cycle_start = 0.5 * p.fs
            offset = floor((episode - 1) * p.fs + 1);
            index1 = offset + cycle_start;
            index2 = index1 + floor(7 * p.fs / p.bias_modulation_rate_C);
            uncorrelated_rate3(index1:index2) = uncorrelated_rate3(index1:index2) .* sine_mod(1:index2-index1+1);
            correlated_rate3(index1-coffset:index2-coffset) = correlated_rate3(index1-coffset:index2-coffset) .* sine_mod(1:index2-index1+1);
        end
    end
end

% Cycle triggered stimulation based on LFPB threshold crossings.
p.stim_phase_sine = sin(p.stim_phase * pi / 180); % Used to calculate amplitude level where stimulation should occur.
if (p.stim_phase < 0)
    p.stim_phase = p.stim_phase + 360;
end
if p.conditioning_type == 3
    if (p.stim_phase < 90) || (p.stim_phase >= 270)
        p.lfp_detect_level = -p.lfp_detect_level; % Trigger on rising through an amplitude level after falling below negative detect level
    end
end

%%%%

Init_Bias(1, p.bias_value, correlated_rate1); % Correlated bias probability for a column is always first.
Init_Bias(2, p.bias_value, uncorrelated_rate1); % Uncorrelated bias probability for column must come immediately after.
Init_Bias(3, p.bias_value, correlated_rate2);
Init_Bias(4, p.bias_value, uncorrelated_rate2);
Init_Bias(5, p.bias_value, correlated_rate3);
Init_Bias(6, p.bias_value, uncorrelated_rate3);
Init_Bias(7, p.bias_value, 0); % No correlated biases for outputs
Init_Bias(8, p.bias_value, out_rate);

% Fill in unit information

Init_Unit('Ae', p.n_excit, p.threshold, 0, 1, 1);  % Column A, Bias index for Col A, LFP index 1
Init_Unit('Ai', p.n_inhib, p.threshold, 0, 1, 1);
Init_Unit('Ao', p.n_out, p.output_threshold_range, 0 + 16384, 7, 11); % Output units tagged wit +16384

Init_Unit('Be', p.n_excit, p.threshold, 1, 3, 2);  % Column B, Bias index for Col B, LFP index 2
Init_Unit('Bi', p.n_inhib, p.threshold, 1, 3, 2);
Init_Unit('Bo', p.n_out, p.output_threshold_range, 1 + 16384, 7, 12);

Init_Unit('Ce', p.n_excit, p.threshold, 2, 5, 3);  % Column C, Bias index for Col C, LFP index 3
Init_Unit('Ci', p.n_inhib, p.threshold, 2, 5, 3);
Init_Unit('Co', p.n_out, p.output_threshold_range, 2 + 16384, 7, 13);

% Fill in connection information

Init_Weight('Ae', 1:p.n_excit, 'Ae', 1:p.n_excit, p.epsp2excit_incol);
Init_Weight('Ae', 1:p.n_excit, 'Ai', 1:p.n_inhib, p.epsp2inhib_incol);
Init_Weight('Ae', 1:p.n_excit, 'Be', 1:p.n_excit, p.epsp2excit_outcol);
Init_Weight('Ae', 1:p.n_excit, 'Bi', 1:p.n_inhib, p.epsp2inhib_outcol);
Init_Weight('Ae', 1:p.n_excit, 'Ce', 1:p.n_excit, p.epsp2excit_outcol);
Init_Weight('Ae', 1:p.n_excit, 'Ci', 1:p.n_inhib, p.epsp2inhib_outcol);
Init_Weight('Ai', 1:p.n_inhib, 'Ae', 1:p.n_excit, p.ipsp2excit_incol);
Init_Weight('Ai', 1:p.n_inhib, 'Ai', 1:p.n_inhib, p.ipsp2inhib_incol);
Init_Weight('Ae', 1:p.n_excit, 'Ao', 1:p.n_out, p.epsp2output_incol);

Init_Weight('Be', 1:p.n_excit, 'Be', 1:p.n_excit, p.epsp2excit_incol);
Init_Weight('Be', 1:p.n_excit, 'Bi', 1:p.n_inhib, p.epsp2inhib_incol);
Init_Weight('Be', 1:p.n_excit, 'Ae', 1:p.n_excit, p.epsp2excit_outcol);
Init_Weight('Be', 1:p.n_excit, 'Ai', 1:p.n_inhib, p.epsp2inhib_outcol);
Init_Weight('Be', 1:p.n_excit, 'Ce', 1:p.n_excit, p.epsp2excit_outcol);
Init_Weight('Be', 1:p.n_excit, 'Ci', 1:p.n_inhib, p.epsp2inhib_outcol);
Init_Weight('Bi', 1:p.n_inhib, 'Be', 1:p.n_excit, p.ipsp2excit_incol);
Init_Weight('Bi', 1:p.n_inhib, 'Bi', 1:p.n_inhib, p.ipsp2inhib_incol);
Init_Weight('Be', 1:p.n_excit, 'Bo', 1:p.n_out, p.epsp2output_incol);

Init_Weight('Ce', 1:p.n_excit, 'Ce', 1:p.n_excit, p.epsp2excit_incol);
Init_Weight('Ce', 1:p.n_excit, 'Ci', 1:p.n_inhib, p.epsp2inhib_incol);
Init_Weight('Ce', 1:p.n_excit, 'Be', 1:p.n_excit, p.epsp2excit_outcol);
Init_Weight('Ce', 1:p.n_excit, 'Bi', 1:p.n_inhib, p.epsp2inhib_outcol);
Init_Weight('Ce', 1:p.n_excit, 'Ae', 1:p.n_excit, p.epsp2excit_outcol);
Init_Weight('Ce', 1:p.n_excit, 'Ai', 1:p.n_inhib, p.epsp2inhib_outcol);
Init_Weight('Ci', 1:p.n_inhib, 'Ce', 1:p.n_excit, p.ipsp2excit_incol);
Init_Weight('Ci', 1:p.n_inhib, 'Ci', 1:p.n_inhib, p.ipsp2inhib_incol);
Init_Weight('Ce', 1:p.n_excit, 'Co', 1:p.n_out, p.epsp2output_incol);

% Finish initialzation. Sort weights for speed, etc.

Init_Done();

% Setup output directory

ct = datevec(now);
fileindex = 1;
datestamp = [num2str(ct(1)) num2str(ct(2),'%02d') num2str(ct(3),'%02d') '_' num2str(fileindex,'%02d')];
foldername = sprintf('%s\\%s\\%s_%s', p.folder, p.subfolder, p.prefix, datestamp);
while exist(foldername, 'file')
    fileindex = fileindex + 1;
    datestamp = [num2str(ct(1)) num2str(ct(2),'%02d') num2str(ct(3),'%02d') '_' num2str(fileindex,'%02d')];
    foldername = sprintf('%s\\%s\\%s_%s', p.folder, p.subfolder, p.prefix, datestamp);
end

if ~exist(foldername, 'dir')
    [status,message,messageid] = mkdir(foldername); %#ok<ASGLU>
    if status == 0
        warndlg(['Could not create the data folder ' foldername],'Save Error');
        return;
    end
else
    warndlg(['Cannot overwrite an existing data folder' foldername],'Save Error');
    return;
end 

disp(['Saving output files to ' foldername]);
copyfile([p.scriptname '.m'], [foldername '\' p.scriptname '.m'], 'f');
copyfile([p.mexname '.c'], [foldername '\' p.mexname '.c'], 'f');
copyfile([p.mexname '.mexw64'], [foldername '\' p.mexname '.mexw64'], 'f');

% Setup space for spike triggered averages of LFP
% LFP(1..3,:) is excitatory LFP contribution for Col A .. Col C
% LFP(4..6,:) is inhibitory LFP contribution for Col A .. Col C
% LFP(7,:) is conditioning type dependent.  For Gamma triggered stimulation, it is the gamma filtered LFP A. For EMG triggered it is the band pass filtered EMG A
% LFP(8..10,:) is synthethic EMG A .. EMG C
% LFP(11..13,:) is output LFP for each column. (not currently used for anything)

sweeps = zeros(p.time_steps, p.units); % spike histogram for each unit.
sweep_count = 0;
activity = zeros(p.units, p.time_steps, 'double'); % Save 10 seconds of activity for each unit.
p.lfp_count = max(p.unit_lfp_offset) + 1;  % If there are inhibitory connections to output units, then use + 4 (since iLFP is stored 3 indexes beyond eLFP)
lfp = zeros(p.lfp_count, p.time_steps, 'double');  % LFP separated by unit type and column.
p.lfp_cycle_tests = zeros(0, 3);
evoked_potentials = zeros(9, 0);  % Evoked potential (max of test range LFP - mean of baseline LFP).

% This is our progress figure.
hfig = figure;
set(hfig, 'Position', [100 100 900 1000]);

sta_sweeps = zeros(p.units,1);
scalems = [-50 100]; % 50 ms before trigger, 100 ms after.  Use integers.
nbins = (scalems(2) - scalems(1)) * p.msfs;
prebins = -scalems(1) * p.msfs;
sta = zeros(p.units, nbins);     % spike triggered averages
scalex = (-prebins:nbins-prebins-1) / p.msfs;   % millisecond x-axis scale for averages

% Peri-event time histograms
psthdefs = { ...
    'Ae -> Ae'; 'Ae -> Ai'; 'Ae -> Be'; 'Ae -> Bi'; 'Ae -> Ce'; 'Ae -> Ci'; 'Ai -> Ae'; 'Ai -> Ai'; ...
    'Be -> Ae'; 'Be -> Ai'; 'Be -> Be'; 'Be -> Bi'; 'Be -> Ce'; 'Be -> Ci'; 'Bi -> Be'; 'Bi -> Bi'; ...
    'Ce -> Ae'; 'Ce -> Ai'; 'Ce -> Be'; 'Ce -> Bi'; 'Ce -> Ce'; 'Ce -> Ci'; 'Ci -> Ce'; 'Ci -> Ci'; ...
    ' T -> Ae'; ' T -> Ai'; ' T -> Be'; ' T -> Bi'; ' T -> Ce'; ' T -> Ci'; ' T -> Ao'; ' T -> Bo'; ' T -> Co'}; 
npsth = length(psthdefs);
psthscale = -100:0.5:100; % Scale in milliseconds
psthedges = psthscale * p.msfs; % Bin edges in samples
psth = zeros(length(psthscale), npsth+10); % 10 extra PSTH for ST Ae1->Be and Ae1->Bi, or GT -> Ae/Ai/Be/Bi/Ce/Ci
psthn = zeros(1, npsth+10);                % Number of sweeps in each PSTH (used for converting counts to Hz)
weightedges = (-p.max_strength:10:p.max_strength) + 5;
weightedgespos = (0:10:p.max_strength) + 5;

% List of all weights connecting Column A -> Column B
p.weights_A_to_B = uint32(find((p.unit_column_id(p.weight_pre_unit) == 0) & (p.unit_column_id(p.weight_post_unit) == 1)));
p.sumWeights_A_to_B = zeros(1,p.time_steps); % Simulation sums weights A->B for specific training types.

% Clear return values, setup for first time step.

sequence_start_nocond = 1;
sequence_start_test1 = sequence_start_nocond + p.ntrials_off;
sequence_start_cond = sequence_start_test1 + p.ntrials_test;
sequence_start_test2 = sequence_start_cond + p.ntrials_on;
sequence_end = sequence_start_test2 + p.ntrials_test;
sequence_length = p.ntrials_off + p.ntrials_test + p.ntrials_on + p.ntrials_test;
p.ntrials = p.ntrials_repeat * sequence_length;  % Total number of trials to do
p.sequence_counter = 0;
p.repeat_counter = 0;
for itrial = 1:p.ntrials
    p.itrial = itrial;
    
    p.last_section_trial = 0;  % Flag to indicate last trial in the current sequence section.
    iseq = mod(p.itrial - 1, sequence_length) + 1;
    if iseq == sequence_start_nocond
        p.repeat_counter = p.repeat_counter + 1;
        p.conditioning_flag = 0;
        p.train_on = 1;
        p.sequence_counter = 1;
        section_label = 'Training no conditioning';
    elseif iseq == sequence_start_test1
        p.conditioning_flag = 0;
        p.train_on = 0;
        p.sequence_counter = 1;
        section_label = 'Testing before conditioning';
    elseif iseq == sequence_start_cond
        p.conditioning_flag = 1;
        p.train_on = 1;
        p.sequence_counter = 1;
        section_label = 'Training with conditioning';
        % Save test stim triggered LFP from previous section
        avetestlfp = avelfp;
        avetestlfp_sweeps = avelfp_sweeps;
        statest = sta;
        statest_sweeps = sta_sweeps;
    elseif iseq == sequence_start_test2
        p.conditioning_flag = 0;
        p.train_on = 0;
        p.sequence_counter = 1;
        section_label = 'Testing after conditioning';
    else
        p.sequence_counter = p.sequence_counter + 1;
        nseq = iseq + 1; % Next trial sequence number
        if (nseq == sequence_start_test1) || (nseq == sequence_start_cond) || (nseq == sequence_start_test2) || (nseq == sequence_end)
            p.last_section_trial = 1;
        end
    end
    
    if p.sequence_counter == 1  % Restart averages at the begining of each sequence section.
        sta_sweeps = zeros(p.units,1);
        sta = zeros(p.units, nbins); 
        avelfp = zeros(40,nbins);
        avelfp_sweeps = 0;
        psth = zeros(length(psthscale), npsth+10);
        psthn = zeros(1, npsth+10);
        sweeps = zeros(p.time_steps, p.units);
        sweep_count = 0;
        lfpstore = zeros(p.n_cols, 0);
        p.sumWeights_A_to_B = zeros(1,p.time_steps); % Simulation sums weights A->B for specific training types.
    end

    if (p.conditioning_type == 1)      
        p.rec_unit = p.trigger_unit; % Spike triggered conditioning is on
    else
        p.rec_unit = 0; % No spike trigger unit in other conditioning method
    end
    
    if (p.conditioning_type == 4) || (p.conditioning_type == 5) || (p.conditioning_type == 6)
        % Tetanic conditioning is a fixed chance at each timestep
        chance = p.tetanic_freq / p.fs;
        tstamps = find(random('uniform', 0, 1, [1 floor(p.fs * p.conditioning_secs)]) < chance);
        short = find(diff(tstamps) < p.stim_refractory);
        if ~isempty(short)
            tstamps(short+1) = []; % Remove intervals that are shorter than our refractory period
        end
        p.stim_source_times(:) = 0;
        p.stim_source_times(tstamps) = 1;
        if p.stim_pulse_train >= 2
            % Pulse trains are at 300 Hz.
            for ipulse = 1:p.stim_pulse_train-1
                p.stim_source_times(tstamps + floor(3.3 * p.msfs * ipulse)) = 1;
            end
        end
    end

    % Run network for 1 iteration. Remember clock ticks of any conditioning
    % stimulations for later.
    p.stim_output_times = zeros(p.time_steps, 1, 'uint16'); % Will record the timesteps where conditioning stimuli were delivered.
    p.trigger_times = zeros(p.time_steps, 1, 'uint16');     % Clockticks where the conditioning triggers were detected (even if conditioning is turned off)
    [iUnit, tSpike] = spikenet50mex(p, activity, lfp);
    trig_timestamps = find(p.trigger_times == 1);           % Convert to list of timestamps   
    short = find(diff(trig_timestamps) < p.stim_refractory);%   .. and remove short intervals (stim train intervals are 3.3 ms)
    if ~isempty(short)
        trig_timestamps(short+1) = [];
    end
    stim_timestamps = find(p.stim_output_times == 1);
    clockticks = int32(stim_timestamps + (p.itrial - 1) * p.time_steps); 
    p.stim_output_clockticks = [p.stim_output_clockticks; clockticks];
    for iu = 1:p.units
        % Keep track of number of output spikes for each unit
        p.unit_spike_counts(iu) = p.unit_spike_counts(iu) + length(find(iUnit == iu));
    end
    
    % Calculate PSTHs
    % Trigger PSTHs are always calculated.  Group PSTHs are only calculated
    % for the final cycle (may need to extend this to several of the last
    % cycles for smoother graphs)
    for ipsth = 1:npsth
        psth_name = psthdefs{ipsth}; % Text name of psth
        targ_name = psth_name(end-1:end);  % Target name (eg 'Ae')
        ref_name = psth_name(1:2);   % Reference name (eg 'Ae')

        % Get spike times for reference unit and target units
        if ref_name(2) == 'T'  % Special case for conditioning trigger times
            ref_times = trig_timestamps;
        else
            if (iseq >= sequence_start_test2) && (p.ntrials_repeat == p.repeat_counter)   % && (p.last_section_trial == 1)
                ref_times = tSpike(ismember(iUnit, find(strncmp(p.unit_names, ref_name, 2))));
            else
                continue;
            end
        end
        ref_times = ref_times(ref_times < 7 * p.fs);
        targ_index = find(strncmp(p.unit_names, targ_name, 2));
        targ_times = tSpike(ismember(iUnit, targ_index));

        % Sum a sweep for each reference spike time
        nref = length(ref_times);
        psthn(ipsth) = psthn(ipsth) + nref * length(targ_index);
        for iref = 1:nref
            sweep_times = targ_times - ref_times(iref);
            sweep_times = sweep_times((sweep_times >= psthedges(1)) & (sweep_times <= psthedges(end)));
            if ~isempty(sweep_times)
                psth(:,ipsth) = psth(:,ipsth) + histc(sweep_times, psthedges, 1);
            end
        end
    end

    % Plot output activity of representative units.

    for iun = 1:p.units
        spikes = tSpike(iUnit == iun);
        if ~isempty(spikes)
            sweeps(spikes, iun) = sweeps(spikes, iun) + 1; % Works because spikes should never overlap on the same index.
        end
    end
    
    plotnames = {'Ae'; 'Be'; 'Ce'; 'Ao'; 'Bo'; 'Co'};
    lfpnames = {'A'; 'B'; 'C'; 'A-'; 'B-'; 'C-'; 'A0'; 'B0'; 'C0'};
    condnames = {'No Conditioning'; 'ST'; 'PP'; 'CT'; 'Tetanic-A'; 'Tetanic-B'; 'Tetanic-C'; 'ET'; 'GT'; 'AT'};
    sweep_count = sweep_count + 1;
    for iplot = 1:6
        iun = find(strncmp(p.unit_names, plotnames{iplot}, 2));
        if ~isempty(iun)
            subplot(6, 3, iplot * 3 - 2);
            sumswp = sum(sweeps(:, iun), 2) / length(iun);
            sweep = 20 * sum(reshape(sumswp, round(p.fs / 20), round(p.time_steps / (p.fs / 20)))) / sweep_count;
            plot(0:0.05:9.95, sweep);
            xlabel([plotnames{iplot} ' (sec)']);
            ylabel('Hz');
            ylim([0 50]);
            if iplot == 1
                if (p.conditioning_type == 3)
                    titlestr = [p.prefix ' ' datestamp ' (' condnames{p.conditioning_type+1} ' ' num2str(p.stim_phase) ' deg, ' num2str(p.stim_uV) ' uV)'];
                elseif (p.conditioning_type == 8)
                    band = [' ' num2str(p.gamma_band(1)) '-'  num2str(p.gamma_band(2)) ' Hz, '];
                    titlestr = [p.prefix ' ' datestamp ' (' condnames{p.conditioning_type+1} band num2str(p.stim_delay_ms) ' ms)'];
                elseif (p.conditioning_type == 0) || (p.conditioning_type == 4) || (p.conditioning_type == 5) || (p.conditioning_type == 6)
                    titlestr = [p.prefix ' ' datestamp ' (' condnames{p.conditioning_type+1} ')'];
                else
                    titlestr = [p.prefix ' ' datestamp ' (' condnames{p.conditioning_type+1} ' ' num2str(p.stim_delay_ms) ' ms, ' num2str(p.stim_uV) ' uV)'];
                end
                titlestr(titlestr == '_') = '-';
                title(titlestr);
            end
        end
    end

    % Plot spike triggered averages of LFP and a few other plots;
    dosta = 1;
    if dosta > 0
        maxbin = nbins - prebins;
        
        % define which plots to show
        source(1) = find(strncmp(p.unit_names, 'Ae', 2), 1); dest(1) = 1; % Ae1..n -> LFPA
        source(2) = source(1); dest(2) = 2; % LFPB
        source(3) = source(1); dest(3) = 3; % LFPC
               
        source(4) = find(strncmp(p.unit_names, 'Be', 2), 1); dest(4) = 1; % Be1..n -> LFPA
        source(5) = source(4); dest(5) = 2; % LFPB
        source(6) = source(4); dest(6) = 3; % LFPC
        
        source(7) = find(strncmp(p.unit_names, 'Ce', 2), 1); dest(7) = 1; % Ce1..n -> LFPA
        source(8) = source(7); dest(8) = 2; % LFPB
        source(9) = source(7); dest(9) = 3; % LFPC

        % plot spike triggered averages of LFP
        
        [b, a] = butter(1, [10 2500] / (p.fs / 2)); %10Hz to 2500Hz butterworth filter
        for iplot = 1:9
            iu = source(iplot);
            ilfp = dest(iplot);
            
            % Single unit spike triggered averages for first unit in each
            % column.  Leave out seconds 8-10 which are involved in stimulus testing.
            ts1 = tSpike(iUnit == iu);  % Spikes from first unit in source column (e.g Ae1 -> Col B)
            lfppos = filter(b, a, lfp(ilfp, :));
            lfpneg = filter(b, a, lfp(ilfp+3, :));
            for ispk = length(ts1):-1:1
                t1 = ts1(ispk);
                if (t1 > nbins) && (t1 <= 7.5 * p.fs - maxbin)
                    sta_sweeps(iplot) = sta_sweeps(iplot) + 1;
                    sta(iplot, :) = sta(iplot, :) + lfppos((t1 - prebins + 1):(t1 + maxbin));
                    sta(iplot+50, :) = sta(iplot+50, :) + lfpneg((t1 - prebins + 1):(t1 + maxbin));
                end
            end
            
            % Column wide spike triggered averages when training is off.
            if (p.train_on == 0)
                ts1 = tSpike((iUnit >= iu) & (iUnit < iu + p.n_excit)); % Spikes from all excitatory source units beginning with source name
                for ispk = length(ts1):-1:1
                    t1 = ts1(ispk);
                    if (t1 > nbins) && (t1 <= 75000 - maxbin)
                        sta_sweeps(iplot+60) = sta_sweeps(iplot+60) + 1;
                        sta(iplot+60, :) = sta(iplot+60, :) + lfpneg((t1 - prebins + 1):(t1 + maxbin)) + lfppos((t1 - prebins + 1):(t1 + maxbin));
                    end
                end
            end

            % Conditioning triggered averages for LFP Ae,Be,Ce,Ai,Bi,Ci, Aout,Bout,Cout
            % column.  Leave out seconds 8-10 which are involved in stimulus testing.
            ts1 = trig_timestamps((trig_timestamps > nbins) & (trig_timestamps <= 7.5 * p.fs - maxbin));  % Conditioning trigger events
            if iplot <= 6
                lfppos = lfp(iplot, :);
            else
                lfppos = lfp(iplot+3, :);
            end
            for ispk = length(ts1):-1:1
                t1 = ts1(ispk);
                sta_sweeps(iplot+40) = sta_sweeps(iplot+40) + 1;
                sta(iplot+40, :) = sta(iplot+40, :) + lfppos((t1 - prebins + 1):(t1 + maxbin));
            end

            % Select correct plot
            if iplot <= 6
                subplot(6, 3, iplot * 3 - 1);
            else
                subplot(6, 3, (iplot - 6) * 3);
            end
            
            ave = sta(iplot,:) / (1000 * sta_sweeps(iplot)); % mV average
            plot(scalex, ave, 'Color', [0 0 0]);
            
            hold on;
            ave = sta(iplot+50,:) / (1000 * sta_sweeps(iplot));
            plot(scalex, ave, 'Color', [1 0 0]);
            
            ave = (sta(iplot,:) + sta(iplot+50,:)) / (1000 * sta_sweeps(iplot));
            plot(scalex, ave, 'Color', [0 0 1]);
                        
            ylim([-20 20]);
            xlim([scalems(1) scalems(2)]);
            name = p.unit_names{iu};
            xlabel([name ' -> LFP' lfpnames{ilfp} ' sta (ms), n = ' num2str(sta_sweeps(iplot))]);
            ylabel('mV');
            if iplot == 1
                title('STA Unit Spike -> LFP');
            end
            if iplot == 7
                title([section_label ' (' num2str(p.itrial) ')']);
            end
            hold off;
        end

        % Update stimulus triggered LFP averages for test figures
        if (p.train_on == 0)
            avelfp_sweeps = avelfp_sweeps + 1;
            r1 = p.fs / 20; % Number of timesteps equal to 50 milliseconds
            r2 = p.fs / 40; % 25 milliseconds
            iplot = 1;
            if p.last_section_trial
                evoked_potentials(:, end+1) = 0; %#ok<AGROW>
            end
            for icol = 1:3
               lfppos = filter(b, a, lfp(icol, :));
               lfpneg = filter(b, a, lfp(icol+3, :));
               lfpfilt = filter(b, a, lfp(icol, :) + lfp(icol+3, :));
               for itime = 1:3
                    isample = 1 + 8 * p.fs + (itime - 1) * 0.7 * p.fs;  % Test stims are placed at specific times
                    range = round(isample - r1) : round(isample + nbins - r1 - 1);
                    avelfp(iplot+20, :) = avelfp(iplot+20, :) + lfppos(range); % Excite LFP
                    avelfp(iplot+30, :) = avelfp(iplot+30, :) + lfpneg(range); % Inhib LFP
                    avelfp(iplot, :) = avelfp(iplot, :) + lfpfilt(range); % Total LFP
                    avelfp(iplot+10, :) = avelfp(iplot+10, :) + lfp(icol+7, range); % Synthetic EMG from output units                 
                    if p.last_section_trial
                        % At the end of each testing section, calculate evoked
                        % potentials for each column to each other column
                        evoked_potentials(iplot, end) = (max(avelfp(iplot,r1:r1+r2-1)) - mean(avelfp(iplot,r2:r1-1))) / avelfp_sweeps;
                    end
                    iplot = iplot + 1;
                end
            end
             
            % Spike triggered averages for each unit in the Ae group onto LFPB
            lfpfilt = filter(b, a, lfp(2,:) + lfp(2+3,:));
            for iplot = 1:40
                ts1 = tSpike(iUnit == iplot);  % Spikes from unit in source column (e.g Ae1 -> Col B)
                for ispk = length(ts1):-1:1
                    t1 = ts1(ispk);
                    if (t1 > nbins) && (t1 <= 7.5 * p.fs - maxbin)
                        sta_sweeps(iplot+70) = sta_sweeps(iplot+70) + 1;
                        sta(iplot+70, :) = sta(iplot+70, :) + lfpfilt((t1 - prebins + 1):(t1 + maxbin));
                    end
                end
            end         
        end

    end
    
    % Calculate average weight from A -> B
    index = find((p.weight_pre_unit <= p.n_excit) & (p.weight_post_unit > p.n_excit + p.n_inhib + p.n_out) & (p.weight_post_unit <= 2*p.n_excit + 2*p.n_inhib + p.n_out));
    mean_weight_AB = mean(p.weight_strength(index)) / p.psp_factor;
    index = find((p.weight_pre_unit <= p.n_excit) & (p.weight_post_unit > 2*(p.n_excit + p.n_inhib + p.n_out)));
    mean_weight_AC = mean(p.weight_strength(index)) / p.psp_factor;
    
    % Display progress line
    progstr = ['Sweep ' num2str(p.itrial) ' Spikes ' num2str(length(tSpike)) ' Stims ' num2str(p.train_info(3)) ' WeightAB ' num2str(mean_weight_AB) ' WeightAC ' num2str(mean_weight_AC)];
    disp(progstr);
    fid = fopen([foldername '\progress.txt'], 'a');
    fprintf(fid, '%s\r\n', progstr);
    fclose(fid);
    
    if p.conditioning_type == 1
        % For Spike triggered conditioning, plot PSTH of Ae1->Ae
        % Get spike times for reference unit and target units
        ref_times = tSpike(iUnit == p.trigger_unit);
        %ref_times = tSpike(ismember(iUnit, find(strncmp(p.unit_names, %'Ae1', 3), 1))); % Old version only supported Ae1 -> Be
        ref_times = ref_times(ref_times < p.conditioning_secs * 7 * p.fs); % Limit spikes to times before test stims
        targ_times = tSpike(ismember(iUnit, find(strncmp(p.unit_names, 'Ae', 2))));
        ipsth = npsth+1; % The space we allocated for this plot in psth().
        
        % Sum a sweep for each reference spike time
        nref = length(ref_times);
        for iref = 1:nref
            sweep_times = targ_times - ref_times(iref);
            sweep_times = sweep_times((sweep_times >= psthedges(1)) & (sweep_times <= psthedges(end)));
            if ~isempty(sweep_times)
                psth(:,ipsth) = psth(:,ipsth) + histc(sweep_times, psthedges, 1);
            end
        end
        
        subplot(6, 3, 12);
        plot(psthscale,  100 * psth(:, ipsth) / p.sequence_counter / p.n_excit);
        xlim([-50 50]);
        trigname = p.unit_names{p.trigger_unit};
        xlabel(['PSTH ' trigname ' -> Ae']);
        ylabel('Hz');
        
        % PSTH Ae1 -> all Be
        targ_times = tSpike(ismember(iUnit, find(strncmp(p.unit_names, 'Be', 2))));
        ipsth = npsth+2; % The space we allocated for this plot in psth().
        
        % Sum a sweep for each reference spike time
        for iref = 1:nref
            sweep_times = targ_times - ref_times(iref);
            sweep_times = sweep_times((sweep_times >= psthedges(1)) & (sweep_times <= psthedges(end)));
            if ~isempty(sweep_times)
                psth(:,ipsth) = psth(:,ipsth) + histc(sweep_times, psthedges, 1);
            end
        end
        
        subplot(6, 3, 15);
        plot(psthscale,  100 * psth(:, ipsth) / p.sequence_counter / p.n_excit);
        xlim([-50 50]);
        xlabel(['PSTH ' trigname ' -> Be']);
        ylabel('Hz');

    elseif p.conditioning_type == 3   
        r2 = round(0.05 * p.fs); % Graph 50 ms on either side of stim
        r1 = 2 * r2; % Timesteps equal to 100 ms
        sweep_ave = zeros(1,r1+1);
        sweep_index = find(p.stim_output_times == 1);  % Find cycle triggers with full sweeps
        sweep_index = sweep_index((sweep_index > r2) & (sweep_index <  p.time_steps - r2));
        n = length(sweep_index);
        if n > 0
            % Beta band filtered LFPB is stored in lfp(7,:)
            xaxis_ms = (-r2:r2) / p.msfs;
            for isweep = 1:n
                sweep_start = sweep_index(isweep) - r2;
                sweep_ave = sweep_ave + lfp(7, sweep_start:sweep_start+r1);
            end
            subplot(6, 3, 12);
            plot(xaxis_ms, sweep_ave / length(sweep_index));
            xlabel(['Trig -> fLFPB (n = ' num2str(n) ')']);
            
            sweep_index = find(p.stim_target_times == 1); % find test triggers
            sweep_ave = zeros(1,r1+1);
            n = length(sweep_index);
            if n > 0
                for isweep = 1:n
                    sweep_start = sweep_index(isweep) - r2;
                    sweep_ave = sweep_ave + lfp(7, sweep_start:sweep_start+r1);
                end
                subplot(6, 3, 15);
                plot(xaxis_ms, sweep_ave / length(sweep_index));
                xlabel(['Test -> fLFPB (n = ' num2str(n) ')']);
            end
        end
        
    elseif (p.conditioning_type == 7) || (p.conditioning_type == 8)
        % EMG or Gamma triggered LFP-A average and ColA spikes
        r2 = round(0.05 * p.fs); % Graph 50 ms on either side of triggers
        r1 = 2 * r2;  % Timesteps equal to 100 ms
        sweep_ave = zeros(1,r1+1);
        sweep_index = find(p.stim_output_times == 1);  % Find gamma triggers with full sweeps
        sweep_index = sweep_index((sweep_index > r2) & (sweep_index <  p.time_steps - r2));
        n = length(sweep_index);
        if n > 0
            % Average LFPA (this is LFPA_pos + LFPA_neg)
            xaxis_ms = (-r2:r2) / p.msfs;
            for isweep = 1:n
                sweep_start = sweep_index(isweep) - r2;
                sweep_ave = sweep_ave + lfp(1, sweep_start:sweep_start+r1) + lfp(4, sweep_start:sweep_start+r1);
            end
            subplot(6, 3, 12);
            plot(xaxis_ms, sweep_ave / length(sweep_index));
            xlabel(['T -> LFPA (n = ' num2str(n) ')']);
           
            % EMG or Gamma filtered LFPA is returned in lfp(7,:)
            sweep_ave = zeros(1,r1+1);
            for isweep = 1:n
                sweep_start = sweep_index(isweep) - r2;
                sweep_ave = sweep_ave + lfp(7, sweep_start:sweep_start+r1);
            end
            subplot(6, 3, 15);
            plot(xaxis_ms, sweep_ave / length(sweep_index));
            if p.conditioning_type == 7
                xlabel(['T -> EMGA (n = ' num2str(n) ')']);
            else
                xlabel(['T -> GammaA (n = ' num2str(n) ')']);
            end
        end
    end

    % Replace last plot with a histogram of the weight distribution
    subplot(6, 3, 18);
    hist(double(p.weight_strength((p.weight_strength ~= 0) & (p.weight_training_rule ~= 0))) / p.psp_factor, weightedges);
    xlim([-p.max_strength p.max_strength]);
    ylim([0 500]);
    xlabel('Trainable weight distribution');
    drawnow;
    
    % Accumulate lfp for a limited number of iterations
    % This is only for the PSD plot shown below the raster plot
    % Reduce the iterations allowed if space becomes an issue.
    if p.do_raster_figure && p.sequence_counter <= 200
        for ilfp = 1:p.n_cols
            lfpcol = lfp(ilfp,:) + lfp(ilfp+3,:);
            samples = length(lfpcol);
            lfpstore(ilfp,end+1:end+samples) = lfpcol;
        end
    end
    
    % Save averages on last trial of each test section
    if (p.last_section_trial == 1) 
        meanWeightAB = p.sumWeights_A_to_B ./ p.psp_factor ./ p.sequence_counter;
        if (p.train_on == 0)
            if strcmp(section_label, 'Testing before conditioning')
                saveas(gcf, [foldername '\averages_pre_' num2str(p.itrial) '.fig']);
                save([foldername '\meanWeightAB_pre.mat'], 'meanWeightAB');
            else
                saveas(gcf, [foldername '\averages_post_' num2str(p.itrial) '.fig']);
                save([foldername '\meanWeightAB_post.mat'], 'meanWeightAB');
            end
        else
            if p.conditioning_flag
                saveas(gcf, [foldername '\averages_cond_' num2str(p.itrial) '.fig']);
                save([foldername '\meanWeightAB_cond.mat'], 'meanWeightAB');
            else
                saveas(gcf, [foldername '\averages_nocond_' num2str(p.itrial) '.fig']);
                save([foldername '\meanWeightAB_nocond.mat'], 'meanWeightAB');
            end
        end
        
        % Raster plot.  Plot spike times with LFPs and LFP Power spectrum.
        if p.do_raster_figure
            h = figure;
            position = get(h, 'Position');
            set(h, 'Position', position .* [1 1 1.5 2]);
            subplot(2,1,1);
            
            samples = 1:floor(2*p.fs); % Plot first 2 seconds
            h = [0 0 0]; % Fill with handles to graphics appearing in the Legend.
            colors = {'b', 'r', [0.7 0.7 0.5], 'c', 'm', 'y'};
            lfpreserve = p.n_excit;
            ncolunits = p.n_excit + p.n_inhib + p.n_out;
            yreserve = lfpreserve + ncolunits + 20; % space to reserver for each column
            ystart = 0;
            yticks = [];
            for ilfp = 1:p.n_cols
                colunits = (((ilfp - 1) * ncolunits) + 1) : (ilfp * ncolunits);
                sindex = find((tSpike < samples(end)) & ismember(iUnit, colunits));
                iu = iUnit(sindex);
                miniu = min(iu);
                iui = find(iu >= miniu + p.n_excit + p.n_inhib);
                iu(iui) = iu(iui) + 20;
                iui = find(iu > miniu + p.n_excit);
                iu(iui) = iu(iui) + 10;
                iu = iu - colunits(1) + ystart + lfpreserve + 10;
                plot(tSpike(sindex) ./ p.fs, iu, 'k.');
                hold on;
                lfpcol = lfp(ilfp,samples) + lfp(ilfp+3,samples);
                lfpcol = lfpcol - min(lfpcol);
                lfpcol = (p.n_excit * lfpcol / max(lfpcol)) + ystart;
                h(ilfp) = plot((0:length(lfpcol)-1) ./ p.fs, lfpcol, 'Color', colors{ilfp}, 'LineWidth', 2);
                yticks(end+1) = ystart + floor(p.n_excit/2);
                yticks(end+1) = ystart + floor(p.n_excit + p.n_excit/2 + 10);
                yticks(end+1) = ystart + floor(2*p.n_excit + p.n_inhib/2 + 20);
                yticks(end+1) = ystart + floor(2*p.n_excit + p.n_inhib + p.n_out/2 + 40);
                ystart = ystart + yreserve + 40;
            end
            set(gca, 'fontsize', 16);
            set(gca, 'fontweight', 'bold');
            title(['Spikes: ' section_label ' (' num2str(p.itrial) ')']);
            xlabel('Time(Sec)');
            set(gca, 'YTick', yticks);
            set(gca, 'YTickLabel', {'LFPA', 'Ae', 'Ai', 'Am', 'LFPB', 'Be', 'Bi', 'Bm', 'LFPC', 'Ce', 'Ci', 'Cm'});
            % legend(h, {'LFP A', 'LFP B', 'LFP C'});
            xlabel('Time (Sec)');
            ylim([0 max(iu)+10]);
            xlim([0 samples(end)/p.fs]); % Show first 1 second
            
            subplot(2,1,2); %Spectral Density of LFP
            freq = 0:1:199;
            periodogram_type = 'Power'; % 'Power' or 'PSD'
            for ilfp = 1:p.n_cols
                %[P, F] = periodogram(lfpstore(ilfp, :) - mean(lfpstore(ilfp, :)), [], freq, p.fs, periodogram_type);
                [P,F] = pwelch(lfpstore(ilfp, :) - mean(lfpstore(ilfp, :)),[],[],freq,p.fs); % Usign pwelch over periodogram
                plot(freq, P, 'Color', colors{ilfp}, 'LineWidth', 2);
                %P = sum(reshape(P, 4, 50));
                %plot((0:49)*4 + 2, P, 'Color', colors{ilfp}, 'LineWidth', 2);
                hold on;
            end
            set(gca, 'fontsize', 16);
            set(gca, 'fontweight', 'bold');
            title('LFP Power Spectrum');
            ylabel(periodogram_type);
            xlabel('Hz');
            legend({'LFP A', 'LFP B', 'LFP C'});
            saveas(gcf, [foldername '\Raster-' num2str(p.itrial) '.fig']);
            close(gcf);
        end
        
        % Save PSTH figure for all T -> Column Groups
        psthfig = figure;
        set(psthfig, 'Position', [100 100 900 1000]);
        lfpnames = {'LFP A', 'LFP A', 'LFP B', 'LFP B', 'LFP C', 'LFP C', 'LFP AOut', 'LFP BOut', 'LFP COut'};
        for iplot = 1:9
            % Plot PSTH
            subplot(5,2,iplot);
            data = psth(:, iplot + 24) * 1000 / psthn(iplot + 24) / (psthscale(2) - psthscale(1)); % Normalize to spikes per second
            %bar(psthscale + 0.5 * (psthscale(2) - psthscale(1)), data, 1);  % Use this for solid fill bars
            [xscale, yvals] = skyline(psthscale, data);  % Create a skyline plot for the histogram
            plot(xscale, yvals, 'LineWidth', 2);
            ylabel('Hz');
            ylim([0 60]);
            xlim([-50 50]);
            xticks(-50:10:50);
            xticklabels({'-50', '', '', '', '', '0' , '', '', '', '', '50'});
            xlabel([psthdefs{iplot + 24} ', ' lfpnames{iplot}]);
            if iplot == 1
                title(['PSTH ' char(titlestr)]);
            end
            if iplot == 2
                title([section_label ' (' num2str(p.itrial) ')']);
            end
            
            % Plot LFP of target column group
            if iplot <= 6
                ista = floor((iplot+1) / 2) + 40;  % Plot same LFP to both T->Ae and T->Ai graphs
                data = (sta(ista, :) + sta(ista+3, :)) / 1000 / sta_sweeps(ista); % Sum exitatory and inhibitory contributions to the LFP.
            else
                data = sta(iplot+40, :) / 1000 / sta_sweeps(iplot+40);
            end
            data = data - min(data);
            hold on;
            plot(scalex, data, 'Color', [0.5,0.5,0.5], 'LineWidth', 2);

        end
        drawnow;
        saveas(gcf, [foldername '\Cond_Trig_PSTH-' num2str(p.itrial) '.fig']);
        close(gcf);
    end

end

p.stop_time = clock;
disp(['Done. Run Time (secs) = ' num2str(etime(p.stop_time, p.start_time))]);
save([foldername '\param.mat'], 'p');
save([foldername '\output_spikes.mat'], 'iUnit', 'tSpike'); % Save last 10 seconds of spikes
clockticks = p.stim_output_clockticks;  
save([foldername '\stim_output_clockticks.mat'], 'clockticks'); % Save timestep course of all stimulation
%save([foldername '\activity.mat'], 'activity'); % save activity of all units if desired (this is a large file)
%save([foldername '\lfp.mat'], 'lfp');

% Save spike triggered averages of lfp to text files  This includes the
% single unit Ae1->lfpB through Ae40->lfpB
fid = fopen([foldername '\Ae1-40_lfpB_postcond_STA.txt'], 'w');
if fid > -1
    fprintf(fid, 'ms');
    for iplot = 1:40  % Header label for each spread sheet column.
        name = p.unit_names{iplot};
        label = ['STA ' name ' lfp' lfpnames{2} ' (' num2str(sta_sweeps(iplot+70)) ')'];
        fprintf(fid, '\t%s', label);
    end
    for iline = 1:length(scalex)
        fprintf(fid, '\n%g', scalex(iline));
        for iplot = 1:40
            ave = sta(iplot+70,iline) / (1000 * sta_sweeps(iplot+70));
            fprintf(fid, '\t%g', ave);
        end
    end
    fclose(fid);
end
fid = fopen([foldername '\Ae1-40_lfpB_precond_STA.txt'], 'w');
if fid > -1
    fprintf(fid, 'ms');
    for iplot = 1:40  % Header label for each spread sheet column.
        name = p.unit_names{iplot};
        label = ['STA ' name ' lfp' lfpnames{2} ' (' num2str(statest_sweeps(iplot+70)) ')'];
        fprintf(fid, '\t%s', label);
    end
    for iline = 1:length(scalex)
        fprintf(fid, '\n%g', scalex(iline));
        for iplot = 1:40
            ave = statest(iplot+70,iline) / (1000 * statest_sweeps(iplot+70));
            fprintf(fid, '\t%g', ave);
        end
    end
    fclose(fid);
end

%Plot PSTH
psthfig = figure;
set(psthfig, 'Position', [100 100 900 1000]);
iplot = 1;
for icol = 1:3
    for ipsth = 1:8
        subplot(8,3, (ipsth -1) * 3 + icol);
        vals = psth(:, iplot) * 1000 / psthn(iplot) / (psthscale(2) - psthscale(1));
        zerobin = find(psthedges == 0);
        %vals(zerobin) = 0.5 * (vals(zerobin-1) + vals(zerobin+1)); % average over center bin.
        plot(psthscale, vals);
        xlim([-25 25]);
        ylim([0 25]);
        xlabel(psthdefs{iplot});
        if iplot == 1
            titlestr2 = ['/' p.subfolder '/' p.prefix ' PSTH ' datestamp];
            titlestr2(titlestr2 == '_') = '-';
            title(titlestr2);
        end
        iplot = iplot + 1;
        hold on;
    end
end
drawnow
saveas(gcf, [foldername '\PSTH.fig']);

%Plot average lfp around stimulation times.
%Evoked potentials are maximum value (in response range) - base value
% (just before any possible response).  Within column EPs are probably
% meaningless given that stimulus evoked IPSPs cancel with EPSPs in the LFP
letter = 'ABC';
epfig = figure;
set(epfig, 'Position', [100 100 1000 1000]);
iplot = 1;
avelfp = avelfp(:,:) / (1000 * avelfp_sweeps); % Convert to mV.
avetestlfp = avetestlfp(:,:) / (1000 * avetestlfp_sweeps); % Convert to mV.
r1 = prebins + 1; % Stimulation bin
r2 = r1 + floor(0.025 * p.fs); % end of peak measurment range
lfpscale = (-r1:nbins-r1-1) / p.msfs;
clipvals = zeros(1,7); % Will put relevant MPI values on clipboard for easy paste into Excel
clipind = [7 3 5 1 7 6 2 4 7]; % Indexes where MPI values should go, We don't save same-column values (on the diagonal).
for icol = 1:3
    for itime = 1:3
        subplot(3, 3, iplot)
        plot(lfpscale, [avetestlfp(iplot, :)', avelfp(iplot, :)'], 'linewidth', 3);
        hold on;
%        plot(lfpscale), avelfp(iplot+20, :), 'Color', [0 0 0]); % Excitatory LFP in Black
%        plot(lfpscale, avelfp(iplot+30, :), 'Color', [1 0 0]); % Inhibitory LFP in Red
        epval = max(avelfp(iplot, r1 + p.conduction_delay:r2)) - avelfp(iplot, r1 + p.conduction_delay);
        eptest = max(avetestlfp(iplot, r1 + p.conduction_delay:r2)) - avetestlfp(iplot, r1 + p.conduction_delay);
        xlim([-25 25]);
        label = ['Stim Col ' letter(itime) ' -> LFP ' letter(icol)];
        xlabel(label);
        mpi = (epval - eptest) * 100 / eptest;
        title(['LFP EP ' num2str(epval) ', MPI ' num2str(mpi, 4)]);
        %clipvals(clipind(iplot)) = epval;  % Save absolute EP value to clipboard
        clipvals(clipind(iplot)) = mpi;     % Save MPI to clipboard

        label = ['Stim Col ' letter(itime) ' LFP ' letter(icol)];
        fid = fopen([foldername '\' label '.txt'], 'w');
        if (fid > -1)
            fprintf(fid, ['ms\t' label ' (n=' num2str(avelfp_sweeps) ')\tExcite\tInhib\t']);
            fprintf(fid, ['before cond (n=' num2str(avetestlfp_sweeps) ')\tExcite\tInhib\n']);
            for iline = 1:length(lfpscale)
                fprintf(fid, '%g\t%g\t%g\t%g\t', lfpscale(iline), avelfp(iplot, iline), avelfp(iplot+20, iline), avelfp(iplot+30, iline));
                fprintf(fid, '%g\t%g\t%g\n', avetestlfp(iplot, iline), avetestlfp(iplot+20, iline), avetestlfp(iplot+30, iline));
            end
            fclose(fid);
        end

        iplot = iplot + 1;
    end
end
clipvals(7) = p.train_info(3); % Number of stimulations delivered
clipboard('copy', sprintf('%.4g\t', clipvals)); % Copy MPI values to the clipboard.  This is something I paste into Excel often.

drawnow;
saveas(gcf, [foldername '\lfp_stim_response.fig']);

%Plot average lfp around column trigger times.
%These averages show central correlation peaks which affect the trigger
%evoked potential calculations. These EPs are displayed for information
%but probably should not be used in comparisons.
epfig = figure;
set(epfig, 'Position', [100 100 1000 1000]);
stascale = (-prebins:nbins-prebins-1) / p.msfs; % X axis scale in ms
plotsequence = [1 4 7 2 5 8 3 6 9];
r1 = prebins - floor(0.025 * p.fs) + 1; % 25 ms before trigger time
for iplot = 1:9
    iu = source(iplot);
    ilfp = dest(iplot);
    uname = p.unit_names{iu};
    subplot(3, 3, plotsequence(iplot));
    testlfp = 0.001 * statest(iplot + 60, :) ./ statest_sweeps(iplot + 60);
    postlfp = 0.001 * sta(iplot + 60, :) / sta_sweeps(iplot + 60);
    plot(stascale, [testlfp', postlfp'], 'linewidth', 3);
    hold on;
    % Base measurments are made 25 ms before the trigger to account for the
    % central correlation peak.
    epval = max(postlfp(r1+p.conduction_delay:end)) - postlfp(r1+p.conduction_delay);
    eptest = max(testlfp(r1+p.conduction_delay:end)) - testlfp(r1+p.conduction_delay);
    xlim([scalems(1) 50]);
    label = ['STA ' uname(1) 'e -> LFP ' letter(ilfp)];
    xlabel(label);
    title(['LFP EP ' num2str(epval) ', MPI ' num2str((epval - eptest) * 100 / eptest, 4)]);

    label = ['STA ' uname(1) 'e to LFP' letter(ilfp)];
    fid = fopen([foldername '\' label '.txt'], 'w');
    if (fid > -1)
        fprintf(fid, ['ms\t' label ' (n=' num2str(sta_sweeps(iplot + 60)) ')\tbefore cond (n=' num2str(statest_sweeps(iplot + 60)) ')\n']);
        for iline = 1:length(stascale)
            fprintf(fid, '%g\t%g\t%g\n', stascale(iline), postlfp(iline), testlfp(iline));
        end
        fclose(fid);
    end
end
drawnow;
saveas(gcf, [foldername '\lfp_spike_response.fig']);
close(gcf);

%Plot average synthetic EMG around stimulation times.
epfig = figure;
set(epfig, 'Position', [100 100 1000 1000]);
iplot = 1;
for icol = 1:3
    for itime = 1:3
        subplot(3, 3, iplot)
        plot(lfpscale, [avetestlfp(iplot+10, :)', avelfp(iplot+10, :)'], 'linewidth', 3);
        hold on;
        xlim([-10 scalems(2)]);
        label = ['Stim Col ' letter(itime) ' -> EMG ' letter(icol)];
        xlabel(label);
        label = ['Stim Col ' letter(itime) ' to EMG ' letter(icol)];
        fid = fopen([foldername '\' label '.txt'], 'w');
        if (fid > -1)
            fprintf(fid, ['ms\t' label ' (n=' num2str(avelfp_sweeps) ')\t']);
            fprintf(fid, ['before cond (n=' num2str(avetestlfp_sweeps) ')\n']);
            for iline = 1:length(lfpscale)
                fprintf(fid, '%g\t%g\t%g\n', lfpscale(iline), avelfp(iplot, iline), avetestlfp(iplot, iline));
            end
            fclose(fid);
        end
        iplot = iplot + 1;
    end
end
drawnow;
saveas(gcf, [foldername '\emg_spike_response.fig']);
close(gcf);

% Pre and Post Conditioning Evoked Potentials.  This is more interesting
% when p.ntrials_repeat is > 1 and will show evoked potentials for
% multiple repeats through the training cycles.
fid = fopen([foldername '\Mean_Evoked_Potentials.txt'], 'w');
epfig = figure;
set(epfig, 'Position', [100 100 900 1000]);
if fid > -1
    [rows, cols] = size(evoked_potentials); %#ok<ASGLU>
    n = floor(cols / 2);
    fprintf(fid, 'Repeat\t');
    fprintf(fid, 'Pre(A=>A)\tPost(A=>A)\tPre(B=>A)\tPost(B=>A)\tPre(C=>A)\tPost(C=>A)\t');
    fprintf(fid, 'Pre(A=>B)\tPost(A=>B)\tPre(B=>B)\tPost(B=>B)\tPre(C=>B)\tPost(C=>B)\t');
    fprintf(fid, 'Pre(A=>C)\tPost(A=>C)\tPre(B=>C)\tPost(B=>C)\tPre(C=>C)\tPost(C=>C)\r\n');
    for iline = 1:n
        fprintf(fid, '%d', iline);
        for iplot = 1:9
            fprintf(fid, '\t%g\t%g', round(evoked_potentials(iplot, iline * 2 - 1)), round(evoked_potentials(iplot, iline * 2)));
        end
        fprintf(fid, '\r\n');
    end
    fprintf(fid, '\r\nmean');
    for iplot = 1:9
        x1 = evoked_potentials(iplot, 1:2:n*2-1);
        x2 = evoked_potentials(iplot, 2:2:n*2);
        subplot(3,3,iplot);
        plot([x1',x2'] ./ 1000);  % Convert to mV
        fprintf(fid, '\t%g\t%g', mean(x1), mean(x2));
    end
    fprintf(fid, '\r\np__mpi');
    for iplot = 1:9
        x1 = evoked_potentials(iplot, 1:2:n*2-1);
        x2 = evoked_potentials(iplot, 2:2:n*2);
        val1 = mean(x1);
        val2 = mean(x2);
        mpi = 100 * (val2 - val1) / val1;
        [h,pval] = ttest(x1, x2); %#ok<ASGLU>
        fprintf(fid, '\t%g\t%.4g', pval, mpi);
        subplot(3,3,iplot);
        label = [letter(rem(iplot-1,3)+1) ' -> ' letter(floor((iplot+2) / 3)) ' (MPI ' num2str(mpi,3) ', p ' num2str(pval) ')'];
        xlabel(label);
    end
    fprintf(fid, '\r\n');
    fclose(fid);
    subplot(3,3,1);
    title(titlestr);
    subplot(3,3,2);
    title('Pre vs Post Conditioning EP)');
    drawnow;
    saveas(gcf, [foldername '\Mean_Evoked_Potentials.fig']);
    close(gcf);
end

% Display Weight Matrix.  Save one that includes all units
figure;
wvals = zeros(p.units, p.units);
for iw = 1:p.weights
    wvals(p.weight_pre_unit(iw), p.weight_post_unit(iw)) = p.weight_strength(iw);
end
wvals(wvals == 0) = NaN;
pcolor(wvals/ p.psp_factor);
title(titlestr);
caxis([-500 500]);
colormap jet;
shading flat;
drawnow
saveas(gcf, [foldername '\weight_matrix.fig']);
close(gcf);
% Display labeled figure that excludes output units.
display_weight_figure(titlestr);
saveas(gcf, [foldername '\weight_matrix_labeled.fig']);

% Display weight distributions for excitatory weights from and to units
wdfig = figure;
set(wdfig, 'Position', [100 100 900 1000]);
iplot = 1;
colexcit = {'Ae', 'Be', 'Ce'};
colinhib = {'Ai', 'Bi', 'Ci'};
colall = {'Aei', 'Bei', 'Cei'};
for from_col = 3:-1:1
    for to_col = 1:3
        subplot(3,3,iplot);
        preunits = find(strncmp(p.unit_names, colexcit{from_col}, 2)); % From source column excitatory units
        postunits = find(strncmp(p.unit_names, colexcit{to_col}, 2));  % To dest column excitatory and inhibitory units
        postunits = [postunits find(strncmp(p.unit_names, colinhib{to_col}, 2))];
        values = p.weight_strength(ismember(p.weight_pre_unit, preunits) & ismember(p.weight_post_unit, postunits)) / p.psp_factor;
        q2 = mean(values);
        q1 = mean(values(values < q2));
        q3 = mean(values(values > q2));
        hist(values, weightedgespos);
        maxy = 50;
        ylim([0 maxy]);
        line([q1 q1 nan q2 q2 nan q3 q3 nan], [0 maxy nan 0 maxy nan 0 maxy nan]);
        xlim([0 p.max_strength]);
        xlabel([colexcit{from_col} ' -> ' colall{to_col} ' ' num2str(round(q1)) ' ' num2str(round(q2)) ' ' num2str(round(q3))]);
        if (iplot == 1), title('Col -> Col Weight Distributions'); end
        if (iplot == 6), title(titlestr); end
        iplot = iplot + 1;
    end
end
drawnow;
saveas(gcf, [foldername '\weight_distribution.fig']);

if p.conditioning_type == 7
    % display firing rates for all output units.
    Output_Hz = p.unit_spike_counts(strncmp(p.unit_names, 'Ao', 2)) / ( p.ntrials * (p.time_steps / p.fs));
    figure;
    bar(Output_Hz);
    xlim([1 length(Output_Hz)]);
    xlabel('Output unit firing rates');
    ylabel('Hz');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Init_Network()
        % Start network initialization
        %disp(['Estimated units: ' num2str(p.units) ', weights: ' num2str(p.weights)]);
        p.units = 0;
        p.weights = 0;
        p.biases = 1;
        % Caclulate psp_factor for converting uV into weights for our
        % slow and fast exponential decay functions. This is based on the
        % peak value in the PSP curve calculated with the 0.1 ms timestep.
        % This keeps the areas under PSP curves with smaller timesteps more
        % comparable.
        slow_decay = 1 - 0.1 / p.PSP_slow_time_constant_ms; 
        fast_decay = 1 - 0.1 / p.PSP_fast_time_constant_ms;   
        timestep_index = (0:1000); % Assuming psp maximal value is within 100ms of start. 
        yslow = slow_decay .^ timestep_index;
        yfast = fast_decay .^ timestep_index;
        p.psp_factor = 1 / max(yslow - yfast);  % Multiplier for psp amplitudes to account for exponential decay function.
        p.max_psp_value = p.max_strength * p.psp_factor;
    end

    function Init_Done()
        % Finish network initialization
        % Resize weight arrays to maximum defined weights
        p.weight_pre_unit = p.weight_pre_unit(1:p.weights);
        p.weight_post_unit = p.weight_post_unit(1:p.weights);
        p.weight_strength = p.weight_strength(1:p.weights);
        p.weight_training_rule = p.weight_training_rule(1:p.weights);
        
        % extra arrays needed for acceleration and information.
        p.weight_pre_sort = zeros(p.weights, 1, 'uint32');   % Offsets to weights sorted by pre_unit
        p.weight_post_sort = zeros(p.weights, 1, 'uint32');  % Offsets to weights sorted by post_unit
         
        % Sort weights by pre and post units and figure offsets into sort lists.
        [sorted, indices] = sort(p.weight_pre_unit);
        p.weight_pre_sort = uint32(indices - 1);
        for iux = 1:p.units
            index = find(sorted == iux, 1);
            if ~isempty(index)
                p.unit_pre_offset(iux) = index - 1;
            end
        end

        [sorted, indices] = sort(p.weight_post_unit);
        p.weight_post_sort = uint32(indices - 1);
        for iux = 1:p.units
            index = find(sorted == iux, 1);
            if ~isempty(index)
                p.unit_post_offset(iux) = index - 1;
            end
        end

        p.weight_strength2 = p.weight_strength; % weight strength at next time step.
        
        % Initialize unit flags for paired pulse and tetanic.
        p.unit_stim_source = Set_Unit_Flags(p.unit_stim_source, p.stim_source_def);
        p.unit_stim_target = Set_Unit_Flags(p.unit_stim_target, p.stim_target_def);
        
        disp(['Final units: ' num2str(p.units) ', weights: ' num2str(p.weights)]);
    end

    function psp_value = uv2psp(uV_value)
        % Convert psp microvolt value to our exponential offset value.
        % Randomize if uV_Value contains 2 values.
        if length(uV_value) == 1
            psp_value = round(p.psp_factor * uV_value);
        else
            psp_value = round(p.psp_factor * random('Uniform', uV_value(1), uV_value(2)));
        end
    end

    function Init_Bias(bias_index, value, rate)
        p.bias_strength(bias_index) = uv2psp(value);
        p.bias_chance(bias_index, 1:p.time_steps) = uint32(floor((2^32 - 1) * (rate / p.fs)));
        p.biases = max(bias_index, p.biases);
    end

    function Init_Unit(prefix, count, thresh, col_index, bias_index, lfp_index) % ('Ae', p.An, threshold, column, bias_id, lfp_id);
       for iux = 1:count
            p.units = p.units + 1;
            p.unit_names{p.units} = [prefix num2str(iux)];  % Text name of unit
            p.unit_group_index(p.units) = iux;   % Unit Index number within its group
            if length(thresh) <= 1
                p.unit_threshold(p.units) = thresh(1);    % Each unit has its own threshold
                p.unit_output_field_strength(p.units) = p.output_field_range(1);
            else
                p.unit_threshold(p.units) = thresh(1) + (thresh(2) - thresh(1)) * (iux - 1) / (count - 1); % Ramp between min and max threshold.
                p.unit_output_field_strength(p.units) = p.output_field_range(1) + (p.output_field_range(2) - p.output_field_range(1)) * (iux - 1) / (count - 1); % Ramp between min and max connection strength.
            end
            p.unit_column_id(p.units) = col_index;     % Column index, 0 = Column A, 1 = Column B, 2 = Column C
            p.unit_bias_offset(p.units) = bias_index - 1; % Zero based bias id   
            p.unit_lfp_offset(p.units) = lfp_index - 1;   % zero based lfp id
            p.unit_stim_source(p.units) = 0;
            p.unit_stim_target(p.units) = 0;
        end
    end

    function flagsout = Set_Unit_Flags(flags, def)
        % Converts a cell array of Name;ID_List pairs to a list of unit flags
        for idef = 1:2:length(def)
            name = def{idef};
            index = def{idef+1};
            flags(ismember(p.unit_group_index, index) & strncmp(p.unit_names, name, length(name))') = 1;          
        end
        flagsout = flags;
    end

    function Init_Weight(pre_name, pre_index, post_name, post_index, initial_psp)
        pre_units = find(strncmp(p.unit_names, pre_name, length(pre_name)));
        post_units = find(strncmp(p.unit_names, post_name, length(post_name)));
        npre = length(pre_index);
        npost = length(post_index);
        psp_prob = abs(initial_psp(3));
        psp_value = p.psp_factor * random('Uniform', initial_psp(1), initial_psp(2), npre, npost);
        psp_train = initial_psp(4);
        if initial_psp(3) >= 0
            % Random chance for each connection
            psp_chance = (random('Uniform', 0, 1, npre, npost) <= psp_prob);  % Chance for each connection to exist.
        else
            % Connection from first fraction of the Pre units.
            psp_chance = zeros(npre, npost);
            psp_chance(1:ceil(psp_prob * npre), :) = 1;
        end
        
        if p.test_No_ColA_ColB_connections
            % Leseion A<->B connections if requested.
            % It's important that all random number selection occurs before
            % this so that lesioned networks have comparable initial
            % connections elsewhere.
            if (pre_name(1) == 'A') && (post_name(1) == 'B')
                return; % Don't make A->B connections
            end
            if (pre_name(1) == 'B') && (post_name(1) == 'A')
                return; % Don't make B->A connections
            end
        end

        for g1 = 1:npre
            for g2 = 1:npost
                id1 = pre_index(g1);
                id2 = post_index(g2);
                if ((id1 ~= id2) || ~strcmp(pre_name, post_name)) && psp_chance(g1,g2) % Disallow connections with the same group ID (no self connections and allow for an unconnected testing pair between groups)
                    p.weights = p.weights + 1;
                    pre = pre_units(id1);    % Index of presynaptic unit
                    post = post_units(id2);  % Index of postsynaptic unit
                    p.weight_pre_unit(p.weights) = pre;
                    p.unit_pre_count(pre) = p.unit_pre_count(pre) + 1;
                    p.weight_post_unit(p.weights) = post;
                    p.unit_post_count(post) = p.unit_post_count(post) + 1;
                    p.weight_strength(p.weights) = psp_value(g1,g2);   % Current connection strength
                    p.weight_training_rule(p.weights) = psp_train;     % Learning rule type (0 = none);
                    % Check for special testing conditions.
                    if p.test_A1_lesion
                        % Testing is done with weights from A1 set to zero.
                        if strcmp(pre_name, 'Ae') && (g1 == 1)
                            p.weight_test_lesion(p.weights) = 1;
                        end
                    end
                    if p.test_AB_lesion
                        % Testing is done with A->B wieghts set to zero.
                        if (pre_name(1) == 'A') && (post_name(1) == 'B')
                            p.weight_test_lesion(p.weights) = 1;
                        end
                    end
                    if p.test_AC_lesion
                       % Testing is done with A->C wieghts set to zero.
                       if (pre_name(1) == 'A') && (post_name(1) == 'C')
                            p.weight_test_lesion(p.weights) = 1;
                        end
                    end
                end
            end
        end
    end
 
    function psp = calcpsp(weight, delay, bins, zerobin)
        % Floating point version of target PSP shape
        psp = zeros(bins, 1);
        s = uv2psp(weight);
        f = s;
        for i=floor(zerobin + delay):bins
            psp(i) = s - f;
            s = s * p.psp_slow_decay;
            f = f * p.psp_fast_decay;
        end
    end

    function ysubramp = rampsub(y)
        % subtracts from x a linear fit of x.
        xindex = 1:length(y);
        poly = polyfit(xindex, y, 1);
        ysubramp = y - polyval(poly, xindex);
    end

    function [xscale, yvals] = skyline(xdata, ydata)
        % Creates a skyline plot from single line plot data
        n = 2 * length(ydata);
        yvals = zeros(1, n);
        xscale = zeros(1, n);
        yvals(1:2:n-1) = ydata;
        yvals(2:2:n) = ydata;
        xscale(1:2:n-1) = xdata;
        xscale(2:2:n-2) = xdata(2:end);
        xscale(end) = 2 * xdata(end) - xdata(end-1);
    end

    function display_weight_figure(usetitle)
        % Creates a labeled weight matrix figure for our standard
        % 40 exitatory/inhibory/output unit networks
        wdfig = figure;
        set(wdfig, 'Position', [100 100 590 500]);
        
        nColUnits = p.n_excit + p.n_inhib + p.n_out;  % number of units in a column
        nDispUnits = p.n_excit + p.n_inhib;  % number of units in col to display
        nDisplay = 3 * (p.n_excit + p.n_inhib);  % Ingnore column output units.
        
        a1Index = 1;
        aOutIndex = a1Index + (p.n_excit + p.n_inhib);
        b1Index = nColUnits + a1Index;
        bOutIndex = b1Index + (p.n_excit + p.n_inhib);
        c1Index = nColUnits + b1Index;
        cOutIndex = c1Index + (p.n_excit + p.n_inhib);
        
        wvals = zeros(nDisplay, nDisplay);
        for iw = 1:p.weights
            pre_unit = p.weight_pre_unit(iw);
            post_unit = p.weight_post_unit(iw);
            
            if pre_unit < b1Index
                ix = pre_unit;
                if pre_unit >= aOutIndex
                    ix = 0;
                end
            elseif pre_unit < c1Index
                ix = pre_unit - p.n_out;
                if pre_unit >= bOutIndex
                    ix = 0;
                end
            else
                ix = pre_unit - 2 * p.n_out;
                if pre_unit >= cOutIndex
                    ix = 0;
                end
            end
            
            if post_unit < b1Index
                iy = post_unit;
                if post_unit >= aOutIndex
                    iy = 0;
                end
            elseif post_unit < c1Index
                iy = post_unit - p.n_out;
                if post_unit >= bOutIndex
                    iy = 0;
                end
            else
                iy = post_unit - 2 * p.n_out;
                if post_unit >= cOutIndex
                    iy = 0;
                end
            end
            
            if (ix > 0) && (iy > 0)
                wvals(ix, iy) = p.weight_strength(iw);
            end
        end
        wvals(wvals == 0) = NaN;
        wvals(1, end) = p.max_strength *  p.psp_factor;
        wvals(end, 1) = -p.max_strength *  p.psp_factor;
        h = pcolor(wvals/ p.psp_factor);
        
        for itick = 41:40:241
            line([1 241], [itick itick], 'Color', 'k');
            line([itick itick], [1 241], 'Color', 'k');
        end
        
        set(gca, 'XTick', 21:40:221);
        set(gca, 'YTick', 21:40:221);
        set(gca, 'YTickLabel', {'Ae', 'Ai', 'Be', 'Bi', 'Ce', 'Ci'});
        set(gca, 'XTickLabel', {'Ae', 'Ai', 'Be', 'Bi', 'Ce', 'Ci'});
        
        ylabel('Pre-synaptic Units');
        xlabel('Post-synaptic Units');
        if nargin > 0
            title(usetitle)
        else
            title('Weight Matrix');
        end
        colormap jet;
        shading flat;
        caxis([-500 500]);
        colorbar;
        drawnow
    end

end