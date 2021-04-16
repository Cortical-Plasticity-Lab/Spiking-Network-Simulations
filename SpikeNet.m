classdef SpikeNet < handle
    %SPIKENET Object-oriented implementation of spikenet50 (Larry Shupe)
    %
    % Syntax:
    %   p = SpikeNet(); % Initialize with default parameters
    %
    %   ----------- Original Description ------------------
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
    %   -----------        End           ------------------
    %
    % Object-oriented approach implemented by Max Murphy 4/14/2021
    
    properties
        scriptpath
        scriptname
        % Output is saved to a folder named [folder '\' subfolder '\' prefix '_' datestring '_' indexstring]
        folder                       % Data directory prefix
        subfolder                    % Folder name to hold ouput folders.
        prefix  = 'sn'               % Output folder name prefix.
        mexname = 'spikenet50mex';   % Name of the mex routine that runs the network.
        
        % Simulation parameters.  Note that not all parameters affect all conditioning types.
        conditioning_type = 0;       % 0=No conditioning, 1=spike triggered, 2=paired pulse, 3 = cycle triggered on LFPB, 4 = Tetanic stim on Col A, 5 = Tetanic stim on Col B, 6 = Tetatic stim on Col C, 7 = EMG triggered Stimulation, 8 = Gamma triggered stimulation, 9 = Activity triggered stimulation
        stim_delay_ms = 10;          % Stimulaton delay from rec_unit threshold crossing to start of stimulation PSP in designated stim target units.
        stim_pulse_train = 1;        % Pulses in the 300 Hz stim pulse train for conditioning. Limited to 3 for Paired Pulse stimulation.
        stim_pulse_freq = 300;       % Interstimulus interval for stimlus trains (this will be converted to an integer number of timesteps)
        stim_phase = 0;              % Cycle triggered stim phase, 0 = rising, 90 = peak, 180 = falling, 270 = trough. These phases are fairly exact, others are interpolated from these points. Phases between 46..89 and 316..359 have the most interpolation variability.
        stim_refractory_ms = 10;     % Refractory period on spike triggered and LFP cycle triggered stimulation
        LFP_detect_level = 30000;    % Amplitude of LFP used for cycle triggered conditioning (use 30000 for cycle triggered +/-20%, -20000 or so for gamma triggered depending on %correlated bias drive, 1000 for EMG triggered
        gamma_band = [50 80];        % Bandwidth for conditioning_type in Hz (for example [50 80] for 50Hz to 80Hz bandpass filter.
        EMG_band = [100 2500];       % Column A motor output filter (Hz) for EMG triggered stimulation.
        conditioning_secs = 10;      % Seconds in each trial used for conditioning.  Limits conditioning to beginning of each 10 second trial.
        bias_modulation_rate_A = 0;  % Cycles per second modulation rate on Column A. 0 = no modulation. 19, 20, or 21 for sine wave
        bias_modulation_rate_B = 0;  % Cycles per second modulation rate on Column B. 0 = no modulation. 20 for 20Hz sine wave
        beta_band = [15 25];         % LFP_B filter band (Hz) for cycle triggered stimulation.
        bias_modulation_rate_C = 0;  % Cycles per second modulation rate on Column C. 0 = no modulation. 19, 20 or 21 for sine wave
        bias_modulation_amount = 0.0;  % Size of modulation. 0.2 = +/- 20% of normal rate.  Normally 0. Cycle triggered stim uses 0.2, or set to 0 to just have test pulses for comparison purposes.  Activity dependent stim usues 0.4, and tetanic will use anything as a control for other conditioning types
        bias_modulation_step_secs = 2; % Number of seconds between modulation episodes.  Normally 2, use 5 to test for training decay.
        tetanic_freq = 10;           % Stims per sec for tetanic and activity dependent conditioning, Max stims/sec for exponentially distributed tetanic conditioning (actual rate is always lower because of refractory period)
        test_pulse_train = 1;        % Number of pulses in the stimulus test train
        test_pulse_freq = 300;       % Test pulse rate in Hz (300Hz gives 3.3 ms interval, 500Hz gives 2 ms interval)
        
        test_A1_lesion = 0;          % 1 if testing should be done with weights from unit A1 set to 0.
        test_AB_lesion = 0;          % 1 if testing should be done with weights from Col A to Col B set to 0.
        test_AC_lesion = 0;          % 1 if testing should be done with weights from Col A to Col C set to 0.
        test_No_ColA_ColB_connections = 0; % 0 Normal connectivity, 1 = disallow all connections between Col A and ColB.
        do_raster_figure = 0;        % 1 to calculate LFP power spectrum and show along with a dot raster.  This can be slow.
        
        % File for replaying stimulus events.  Empty string for no replay.
        % This can be used for spike-triggered and cycle-triggered "catch" 
        %   trials.
        % Example:
        %   stim_replay_file = ...
        %    fullfile(pwd, 'spikenet50\sn50_20210415_02', ...
        %                  'stim_output_clockticks.mat');
        stim_replay_file = '';
        
        % Number of trials to run under various conditions.
        ntrials_off = 50;     % Number of initial trials with conditioning off (STDP training on).
        ntrials_on = 50;      % Number of following trials with conditioning on (STDP training on).
        ntrials_test = 50;    % Number of trials for summary figures (both conditioning and STDP training off)
        ntrials_repeat = 1;   % Number of repetitions through the off/test/on/test cycle.
        
        % Network configuration constants.  At this point, these cannot be changed
        % without making changes to most of the analysis and figures.
        n_excit = 40;  % Number of units in each excitatory group
        n_inhib = 40;  % Number of units in each inhibitory grou
        n_out = 40;    % Number of final output units for each Column
        n_cols = 3;    % Number of simulated columns.
        
        % Threshold for cortical and motor output units.
        threshold = 5000;        % uV threshold is currently global for all cortical units
        output_threshold_range = [5000 6000]; % Output units have graded thresholds to simulate connections to muscles.
        output_field_range = [500 1500];      % Output unit graded strengths to the output EMG.
        
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
        
        axonal_delay_ms = 2;     % Axonal mS delay is currently global for all non-output connections.
        dendritic_delay_ms = 1;  % Dendritic mS delay is curretly global for all non-output connections
        output_delay_ms = 10;    % Connection delay for connections to output units
        
        % PSP shape is exp(-t/tau_slow) - exp(-t/tau_fast). This will be converted
        % into a pair of leaky integrators with Euler's method for fast computation
        % using a timestep of 1/fs seconds
        
        fs = 10000;  % Time steps per second. Keep as a multiple of 10000 with a maximum of 100000
        PSP_comparison_window_ms = 0:1000; % Width of window to find peaks for post-synaptic potential comparisons or to apply changes.
        PSP_slow_time_constant_ms = 3.2; % Default synaptic potential shape
        PSP_fast_time_constant_ms = 0.8; % Time constants are in milliseconds
        
        % Bias activity
        biases = 1;
        bias_value = 350;                % uV strenth of bias potientials
        bias_rate = 1800;                % Bias spikes per second to each column unit. Consider smaller values if biases are modulated.
        bias_correlated_fraction = .30;  % Fraction of correlated bias spikes
        correlated_bias_n_std = 8;       % +/- 4 ==> 8 number of standard deviations to limit finite representation of distribution.
        correlated_bias_std_ms = 3;      % <= 5ms. Standard deviation of the normally distributed (correlated) bias spikes (in milliseconds)
        out_bias_rate = 2000;            % Number of uncorrelated bias spikes per second for output units.
        
        % Weight limits.  Weight dependence slows down weight changes as a weight
        % approaches zero or maximum (or -maximum for inhibitory weights).  This is
        % an experimental way to apply some homeostasis to the weights to help keep
        % them more central. Standard simulations use no weight dependence.
        % If weight_dependence is non-zero, strengthening changes are multipled by
        % (1 - (current_weight / weight_limit)) ^ weight_dependence, and weakening
        % changes are multipled by (current_weight / weight_limit) ^ weight_dependence
        
        max_strength = 500;      % Maximum allowed connection strength (uV)
        weight_dependence = 0.0; % 0.001 for almost no weight dependence. Maximum of 1 for linear weight dependence.
        
        % Initial PSP ranges in uV.  Arrays are [uVmin uVmax pChanceForConnection enableTraining]
        initmin = 100;
        initmax = 300;
        EPSP_2_excit_incol = [100 300 1/6 1];     % EPSP for connections to excitatory units within a column
        EPSP_2_inhib_incol = [100 300 1/6 1];     % EPSP for connections to inhibitory units within a column
        EPSP_2_excit_outcol = [100 300 1/6 1];    % EPSP for connections to excitatory units in adjacent columns
        EPSP_2_inhib_outcol = [100 300 1/6 1];    % EPSP for connections to inhibitory units in adjacent columns
        IPSP_2_excit_incol = [-100 -300 1/3 1];   % IPSP for connections to excitatory units within a column
        IPSP_2_inhib_incol = [-100 -300 1/3 1];   % IPSP for connections to inhibitory units within a column
        EPSP_2_output_incol = [350 350 1/3 0];    % EPSP for connections to output units.  These are not usually trained. Negative connection chance means first p * N units rather than random chance.
        
        training_rate = 100;    % Train rate factor for both strengthing (pos) and weakening (neg) sides of the SDTP rule.
        train_weakening = 0.55; % Relative amplitude for the weakening side of the SDTP curve.
        STDP_strengthening_slow_time_constant_ms = 15.4; % STDP strengthening potential shape
        STDP_strengthening_fast_time_constant_ms = 2;
        STDP_weakening_slow_time_constant_ms = 33.3;     % STDP weakening potential shape
        STDP_weakening_fast_time_constant_ms = 2;
        
        % Conditioning related parameters
        stim_uV = 2000      % Amplitude of stimulation pulse (micro-volts)
        test_uV = 3000      % Amplitude of test pulse (micro-volts)
        
        % Reproducibility parameters
        random_seed = 1              % Random generator seed value
        stim_replay_clockticks = -1  % Last value in array must be < 0
    end
    
    properties (Access=public)
        % % These are all set in `initialize_network_parameters` % %
        pair_uV                % Size of the paired-pulse stimulation (second pulse only)
        correlated_bias_rate   % Correlated bias spikes per second delivered to each unit in a column. These will cause common input peaks in cross-correlations of same column units.
        uncorrelated_bias_rate % Number of bias spikes per second for normal firing rate.
        stim_excit_range       % Excitatory units to stimulate
        stim_inhib_range       % Inhibitory units to stimulate
        trigger_unit           % Source unit index for spike triggered stimulation.  0 will disable spike triggered stimulation.
        stim_source_def        % Define stimulation sources for non-delayed stimulation when using
        %                           Tetanic, Cycle-Triggered, Activity-Triggered, and Paired pulse (first
        %                           pulse) conditioning
        stim_target_def        % Define stimulation targets for delay stimulation when using Spike- Cycle-
        %                           EMG- Gamma- Triggered, and Paired Pulse (second pulse) conditioning
        
        % % These are set in `initialize_filter_parameters` % %
        msfs                % Time steps per millisecond.
        timestep2ms         % Milliseconds per timestep
        axonal_delay        % axonal PSP delay in timesteps
        dendritic_delay     % dendridic PSP delay in timesteps
        conduction_delay    % total conduction time in timesteps
        output_delay        % Connection delay for any connection to an output unit.
        stim_delay          % stimulation delay in timesteps
        stim_refractory     % stimulation refractory
        gamma_filter_b      % LFP_A filter transfer function numerator coefficients for gamma triggered stimulation
        gamma_filter_a      % LFP_A filter transfer function denominator coefficients for gamma triggered stimulation
        beta_filter_b       % LFP_B filter transfer function numerator coefficients for cycle triggered stimulation
        beta_filter_a       % LFP_B filter transfer function denominator coefficients for cycle triggered stimulation
        EMG_filter_b        % Motor output filter transfer function numerator coefficients for EMG triggered stimulation
        EMG_filter_a        % Motor output filter transfer function denominator coefficients for EMG triggered stimulation
        
        % % These are set in `initialize_time_factor_parameters` % %
        PSP_slow_decay
        PSP_fast_decay
        train_pos_slow_decay    % Shape of strengthening STDP rule (slow component).
        train_pos_fast_decay    % Shape of strengthening STDP rule (fast component).
        train_neg_slow_decay    % Shape of weakening STDP rule (slow component).
        train_neg_fast_decay    % Shape of weakening STDP rule (fast component).
        train_pos_factor        % STDP Strengthening factor (r) for spike pairs where (tPost - tPre) >= 0.
        train_neg_factor        % STDP weakening factor (cr) for spike pairs where (tPost - tPre < 0).
        
        % % These are set in `initialize_normal_pdf_table` % %
        correlated_bias_max_timesteps
        normal_pdf_table
        c_offset
        
        % % These are set in `initialize_network_units` % %
        units
        weights
        time_steps
        train_info
        bias_strength
        bias_chance
        unit_names
        unit_group_index
        unit_threshold
        unit_bias_offset
        unit_pre_count
        unit_post_count
        unit_pre_offset
        unit_post_offset
        unit_LFP_offset
        unit_column_id
        unit_stim_source
        unit_stim_target
        unit_output_field_strength
        weight_pre_unit
        weight_post_unit
        weight_strength
        weight_training_rule
        weight_test_lesion
        
        % % These are set in `initialize_stim_data` % %
        stim_source_times
        stim_target_times
        stim_test_times
        stim_output_clockticks
        unit_spike_counts
        stim_pulse_isi
        init_test_times
        test_pulse_isi
        
        % % These are set in `init_network_` private function % %
        PSP_factor
        max_PSP_value
        
        % These are set when the simulation is run.
        start_time
        version_id
        random_stream
    end
    
    properties (Hidden, Access=public)
        pair_uV_ = 2000;         % Size of the paired-pulse stimulation (second pulse only)
        stim_excit_range_ = inf  % Can be: {inf (def) => all; [] => none; => array of indices to stim}
        stim_inhib_range_ = inf  % Can be: {inf (def) => all; [] => none; => array of indices to stim}
        weights_ = 50000;        % Initial number of weights, but will be adjusted up or down if necessary.
    end
    
    methods
        function p = SpikeNet(np, fp, tfp, pdfp, weights, init_times, stim_replay_file, varargin)
            %SPIKENET  Constructor method (special).
            %
            % Syntax:
            %   p = SpikeNet();
            %   p = SpikeNet(np, fp, tfp, pdfp, weights, init_times,
            %           stim_replay_file, varargin);
            %
            % Inputs:
            %   np - Network parameters struct or cell array of 'name',
            %           value pairs. See also:
            %               SpikeNet.initialize_network_parameters
            %   fp - Filter parameters struct or cell array of 'name',
            %           value pairs. See also:
            %               SpikeNet.initialize_filter_parameters
            %   tfp - Time-Factor parameters struct or cell array of
            %           'name', value pairs. See also:
            %               SpikeNet.initialize_time_factor_parameters
            %   pdfp - Normal-PDF parameters struct or cell array of
            %           'name', value pairs. See also: 
            %               SpikeNet.initialize_normal_pdf_table
            %   duration - Total duration to run a single iteration of the
            %               network simulation test (seconds). Default
            %               value is 10 seconds.
            %   weights - Initial guess on number of weights (connections)
            %               between units of all columns. See also:
            %               SpikeNet.initialize_network_units
            %   init_times - Start times for test stimulation. See also:
            %               SpikeNet.initialize_stim_data
            %   stim_replay_file - Name of stimulation replay file (default
            %            is ''; only used for "catch" triles on "triggered"
            %            stimulation runs). See also:
            %               SpikeNet.initialize_stim_data
            %
            %   varargin - <'Name', value> property value pairs for any of
            %               the public parameter properties.
            if nargin < 1
               np = p.default_network_parameters();
            end
            if nargin < 2
               fp = p.default_filter_parameters();
            end
            if nargin < 3
               tfp = p.default_time_factor_parameters(); 
            end
            if nargin < 4
               pdfp = p.default_normal_pdf_table_parameters(); 
            end
            if nargin < 5
               duration = 10; 
            end
            if nargin < 6
               weights = 50000; 
            end
            if nargin < 7
               init_times = [8.0, 8.7, 9.4]; 
            end
            if nargin < 8
               stim_replay_file = ''; 
            end
            p = SpikeNet.update_fields_from_args(p, varargin{:});
            if iscell(np)
                p.initialize_network_parameters(np{:});
            else
                p.initialize_network_parameters(np);
            end
            if iscell(fp)
                p.initialize_filter_parameters(fp{:});
            else
                p.initialize_filter_parameters(fp);
            end
            if iscell(tfp)
                p.initialize_time_factor_parameters(tfp{:});
            else
                p.initialize_time_factor_parameters(tfp);
            end
            if iscell(pdfp)
                p.initialize_normal_pdf_table(pdfp{:});
            else
                p.initialize_normal_pdf_table(pdfp); 
            end
            p.initialize_network_units(duration, weights);
            p.initialize_stim_data(init_times, stim_replay_file);
        end
        
        %         function delete(p)
        %             %DELETE Ensure deletion of any associated graphics handles.
        %
        %         end
    end
    
    methods (Access=public)
        function run(p, output_folder)
            %RUN Run the simulation.
            %
            % Syntax:
            %   p.run();
            %   p.run(output_folder);
            %
            % Inputs:
            %   output_folder - (char; optional) Name of output folder
            %       (as a relative filepath relative to the project folder)
            %   -> If it does not exist, then this folder will be created
            %       automatically.
            if nargin < 2
                p.initialize_data_file_locations();
            else
                p.initialize_data_file_locations(output_folder);
            end
            p.init_network_();
            p.start_time = clock();
            p.random_stream = rng(p.random_seed);
            p.version_id = mfilename;
        end
        
        function initialize_data_file_locations(p, output_folder)
            %INITIALIZE_DATA_FILE_LOCATIONS Initialize filename parameters
            %
            % Syntax:
            %   p.initialize_data_file_locations();
            %   p.initialize_data_file_locations(output_folder);
            %
            % Inputs:
            %   output_folder - (char; optional) Name of output folder
            %       (as a relative filepath relative to the project folder)
            %   -> If it does not exist, then this folder will be created
            %       automatically.
            
            NAME_SCHEMA = 'spikenet-%d-%d-%d-%dx';
            SCHEMA_DATA = {p.ntrials_off, ...
                           p.ntrials_on, ...
                           p.ntrials_test, ...
                           p.ntrials_repeat};
            SCHEMA_DESCRIPTORS = {...
                'Number of initial trials', ...
                'Number of conditioning trials', ...
                'Number of test trials (following conditioning)', ...
                'Number of OFF-CONDITIONING-TEST cycle repetitions'};
            to_write = [SCHEMA_DESCRIPTORS; SCHEMA_DATA];
            clc;
            fprintf(1,'<strong>%67s</strong>: <strong>%s</strong>\n', ...
                'Parameter', 'Value');
            fprintf(1,'-> %64s: %d\n', to_write{:});
            
            [p.scriptpath, p.scriptname] = fileparts(mfilename('fullpath'));
            cd(p.scriptpath); % Run inside of this scripts directory.
            if nargin < 2
                p.folder = fullfile(pwd, 'results');
                p.subfolder = fullfile( ...
                    sprintf('Condition-%02d', p.conditioning_type), ...
                    sprintf(NAME_SCHEMA, SCHEMA_DATA{:})...
                );
                output_folder = fullfile(p.folder, p.subfolder);
            else
                [f, sf, ~] = fileparts(output_folder);
                if isempty(f)
                    f = fullfile(pwd, 'results');
                elseif ~contains(f, pwd)
                    f = fullfile(pwd, f);
                end
                if isempty(sf)
                   sf = fullfile( ...
                        sprintf('Condition-%02d', p.conditioning_type), ...
                        sprintf(NAME_SCHEMA, SCHEMA_DATA{:})...
                    );
                end
                p.folder = f;
                p.subfolder = sf;
            end
            if exist(output_folder, 'dir')==0
                mkdir(output_folder);
                fid = fopen(fullfile(output_folder, 'about.txt'), 'w+');
                fprintf(fid,'-> %64s: %d\n', to_write{:});
                fclose(fid);
            end        
        end
        
        function initialize_network_parameters(p,varargin)
            %INITIALIZE_NETWORK_PARAMETERS Set network parameters
            %
            % Syntax:
            %   p.initialize_network_parameters();
            %   p.initilaize_network_parameters('par1name',par1value,...);
            %
            % Inputs:
            %   varargin - (Optional) 'Name', value pairs:
            %               * 'bias_correlated_fraction'
            %               * 'bias_rate'
            %               * 'pair_uV'
            %               * 'stim_excit_range'
            %               * 'stim_inhib_range'
            %
            %   Alternatively, can give `varargin` as a struct that has all
            %   of those fields.
            %
            % See also: SpikeNet.default_network_parameters
            
            % Default values for pre-requisite parameters
            pars = p.default_network_parameters();
            
            % Parse optional (updating) parameters input
            pars = SpikeNet.update_fields_from_args(pars, varargin{:});
            
            p.correlated_bias_rate = pars.bias_correlated_fraction * pars.bias_rate;
            p.uncorrelated_bias_rate = pars.bias_rate - p.correlated_bias_rate;
            
            if p.conditioning_type == 2
                p.pair_uV = pars.pair_uV;
            else
                p.pair_uV = p.stim_uV;
            end
            
            if isinf(pars.stim_excit_range)
                p.stim_excit_range = 1:p.n_excit;
            else
                p.stim_excit_range = pars.stim_excit_range;
            end
            
            if isinf(pars.stim_inhib_range)
                p.stim_inhib_range = 1:p.n_inhib;
            else
                p.stim_inhib_range = pars.stim_inhib_range;
            end
            
            if p.conditioning_type == 1  % Spike Triggered Stimulation
                p.trigger_unit = 1;      % Source unit index for spike triggered stimulation.  0 will disable spike triggered stimulation.
            else
                p.trigger_unit = 0;
            end
            
            % Define stimulation sources for non-delayed stimulation when using
            % Tetanic, Cycle-Triggered, Activity-Triggered, and Paired pulse (first
            % pulse) conditioning
            if (p.conditioning_type == 2) || (p.conditioning_type == 3) || (p.conditioning_type == 4) || (p.conditioning_type == 9)
                p.stim_source_def = {'Ae'; p.stim_excit_range; 'Ai'; p.stim_inhib_range};  % Tetanic A, Paired Pulse, Activity triggered, or Cycle triggered stim
            elseif (p.conditioning_type == 5)
                p.stim_source_def = {'Be'; p.stim_excit_range; 'Bi'; p.stim_inhib_range};  % Tetanic B stim
            elseif (p.conditioning_type == 6)
                p.stim_source_def = {'Ce'; p.stim_excit_range; 'Ci'; p.stim_inhib_range};  % Tetanic C stim
            else
                p.stim_source_def = {}; % No sources to stimulate
            end
            
            % Define stimulation targets for delay stimulation when using Spike- Cycle-
            % EMG- Gamma- Triggered, and Paired Pulse (second pulse) conditioning
            if (p.conditioning_type == 1) || (p.conditioning_type == 2) || (p.conditioning_type == 7) || (p.conditioning_type == 8)
                p.stim_target_def = {'Be'; p.stim_excit_range; 'Bi'; p.stim_inhib_range};
            elseif (p.conditioning_type == 3)
                % Cycle triggered stimulation when using LFPB as trigger and Col A as stimulation target.
                p.stim_target_def = {'Ae'; p.stim_excit_range; 'Ai'; p.stim_inhib_range};
            else
                p.stim_target_def = {}; % No targets to stimulate
            end
        end
        
        function initialize_filter_parameters(p,varargin)
            %INITIALIZE_FILTER_PARAMETERS Set filter pass bands using fs
            %
            % Syntax:
            %   p.initialize_filter_parameters();
            %   p.initilaize_filter_parameters('par1name',par1value,...);
            %
            % Inputs:
            %   varargin - (Optional) 'Name', value pairs:
            %               * 'fs'
            %               * 'axonal_delay_ms'
            %               * 'dendritic_delay_ms'
            %               * 'output_delay_ms'
            %               * 'stim_delay_ms'
            %               * 'stim_refractory_ms'
            %               * 'gamma_band'
            %               * 'beta_band'
            %               * 'EMG_band'
            %
            %   Alternatively, can give `varargin` as a struct that has all
            %   of those fields.
            %
            % See also: SpikeNet.default_filter_parameters
            
            % Default values for pre-requisite parameters
            pars = p.default_filter_parameters();
            
            % Parse optional (updating) parameters input
            pars = SpikeNet.update_fields_from_args(pars, varargin{:});
            
            % Update the relevant fields (using field names of `pars`)
            p = SpikeNet.update_fields_from_struct(p, pars);
            
            % Make sure p.fs is multiple of 10000 between 10000 and 100000.
            p.fs = min(100000, max(10000, round(p.fs / 10000) * 10000));
            p.msfs = round(p.fs / 1000); % Time steps per millisecond.
            p.axonal_delay = floor(p.axonal_delay_ms * p.msfs);        % axonal PSP delay in timesteps
            p.dendritic_delay = floor(p.dendritic_delay_ms * p.msfs);  % dendridic PSP delay in timesteps
            p.conduction_delay = p.axonal_delay + p.dendritic_delay;   % total conduction time in timesteps
            p.output_delay = floor(p.output_delay_ms * p.msfs);        % Connection delay for any connection to an output unit.
            p.stim_delay = floor(p.stim_delay_ms * p.msfs);            % stimulation delay in timesteps
            p.stim_refractory = floor(p.stim_refractory_ms * p.msfs);  % stimulation refractory
            [p.gamma_filter_b, p.gamma_filter_a] = butter(1, p.gamma_band / (p.fs / 2)); % LFP_A filter for gamma triggered stimulation
            [p.beta_filter_b, p.beta_filter_a] = butter(1, p.beta_band / (p.fs / 2)); % LFP_B filter for cycle triggered stimulation
            [p.EMG_filter_b, p.EMG_filter_a] = butter(1, p.EMG_band / (p.fs / 2)); % Motor output filter for EMG triggered stimulation
        end
        
        function initialize_time_factor_parameters(p, varargin)
            %INITIALIZE_TIME_FACTOR_PARAMETERS Set the time decay factors for PSP/STDP
            %
            % Syntax:
            %   p.initialize_time_factor_parameters();
            %   p.initilaize_time_factor_parameters('par1name',par1value,...);
            %
            % Inputs:
            %   varargin - (Optional) 'Name', value pairs:
            %               * 'PSP_slow_time_constant_ms'
            %               * 'PSP_fast_time_constant_ms'
            %               * 'STDP_strengthening_slow_time_constant_ms'
            %               * 'STDP_strengthening_fast_time_constant_ms'
            %               * 'STDP_weakening_slow_time_constant_ms'
            %               * 'STDP_weakening_fast_time_constant_ms'
            %               * 'training_rate'
            %               * 'train_weakening'
            %               * 'train_pos_factor'
            %
            %   Alternatively, can give `varargin` as a struct that has all
            %   of those fields.
            %
            % See also: SpikeNet.default_time_factor_parameters
            
            % Default values for pre-requisite parameters
            pars = p.default_time_factor_parameters();
            
            % Parse optional (updating) parameters input
            pars = SpikeNet.update_fields_from_args(pars, varargin{:});
            
            % Update the relevant fields (using field names of `pars`)
            p = SpikeNet.update_fields_from_struct(p, pars);
            
            % Convert time constants into decay factors for PSP and STDP shapes
            % (this is Eulers method used for estimating an exponential decay function)
            p.timestep2ms = 1000 / p.fs; % Milliseconds per timestep
            p.PSP_slow_decay = 1 - p.timestep2ms / p.PSP_slow_time_constant_ms;
            p.PSP_fast_decay = 1 - p.timestep2ms / p.PSP_fast_time_constant_ms;
            p.train_pos_slow_decay = 1 - p.timestep2ms / p.STDP_strengthening_slow_time_constant_ms;   % Shape of strengthening STDP rule.
            p.train_pos_fast_decay = 1 - p.timestep2ms / p.STDP_strengthening_fast_time_constant_ms;
            p.train_neg_slow_decay = 1 - p.timestep2ms / p.STDP_weakening_slow_time_constant_ms;       % Shape of weakening STDP rule.
            p.train_neg_fast_decay = 1 - p.timestep2ms / p.STDP_weakening_fast_time_constant_ms;
            
            p.train_pos_factor = p.training_rate; % STDP Strengthening factor (r) for spike pairs where (tPost - tPre) >= 0.
            p.train_neg_factor = p.train_weakening * p.train_pos_factor; % STDP weakening factor (cr) for spike pairs where (tPost - tPre < 0).
            
        end
        
        function initialize_normal_pdf_table(p, varargin)
            %INITIALIZE_NORMAL_PDF_TABLE  Initialize the random table
            %
            % Syntax:
            %   p.initialize_normal_pdf_table();
            %   p.initialize_normal_pdf_table('par1name',par1value,...);
            %
            % Inputs:
            %   varargin - (Optional) 'Name', value pairs:
            %               * 'correlated_bias_std_ms'
            %               * 'correlated_bias_n_std'
            %
            %   Alternatively, can give `varargin` as a struct that has all
            %   of those fields.
            %
            % See also: SpikeNet.default_normal_pdf_table_parameters
            
            % Default values for pre-requisite parameters
            pars = p.default_normal_pdf_table_parameters();
            
            % Parse optional (updating) parameters input
            pars = SpikeNet.update_fields_from_args(pars, varargin{:});
            
            % Update the relevant fields (using field names of `pars`)
            p = SpikeNet.update_fields_from_struct(p, pars);
            p.correlated_bias_max_timesteps = p.correlated_bias_n_std * p.correlated_bias_std_ms * p.msfs;
            p.normal_pdf_table = round(random('norm', 0.5 * p.correlated_bias_max_timesteps, p.correlated_bias_std_ms * p.msfs, [100000, 1])); % used for correlated bias spikes
            p.normal_pdf_table((p.normal_pdf_table < 0) | (p.normal_pdf_table > p.correlated_bias_max_timesteps)) = [];  % remove elements more than 4 standard deviations from mean.
            p.normal_pdf_table = uint16(p.normal_pdf_table(1:2^16));
            p.c_offset = round(mean(p.normal_pdf_table));  % clocktick offset for any modulation of correlated bias inputs.  This is just the mean all possible time adjustments.

        end
    
        function initialize_network_units(p, duration, weights)
            %INITIALIZE_NETWORK_UNITS Pre-allocate to save unit data
            %
            % Syntax:
            %   p.initialize_network_units();
            %   p.initialize_network_units(duration, weights);
            %
            % Inputs:
            %   duration - Total duration to run a single iteration of the
            %               network simulation test (seconds). Default
            %               value is 10 seconds.
            %
            %   weights - Default value is 50000, although it is adjusted
            %               as needed by the initialization.
            
            if nargin < 2
               duration = 10; 
            end
            
            if nargin < 3
               weights = 50000; 
            end
            
            % Define the number of biases, units, and weights.
            % Biases are used to give units background activity.  
            % This activity can be modulated by the bias_chance array,
            %   which is assigned to each bias.
            % p.time_steps is the number of time steps to run the network 
            %   on each iteration. This value is usually in the 1 to 10 
            %   second range, but may be set to a time based on the length 
            %   of an experimental trial in a simulated task.
            p.units = p.n_cols * (p.n_excit + p.n_inhib + p.n_out);  % Number of actual units, but will be adjusted up if necessary.
            p.weights = weights;        % Initial number of weights, but will be adjusted up or down if necessary.
            p.time_steps = floor(duration * p.fs); % 10 Seconds of simulation time per iteration.
            p.train_info = zeros(1,4);% Place for keeping track of some training summary results

            % Allocate space for biases.  Biases return a randomized value 
            %   on every access using a 32-bit random number generator R.
            % General algorithm:
            %   If R(access) < bias_chance(step):
            %       return bias_strength 
            %   else 
            %       return 0
            
            % Bias[1] is for uncorrelated biases. 
            % Bias[2, 3, 4] are for correlated biases, assigned to groups 
            %     [A, B, C] respectively.
            % Bias[5] is for output units.
            p.bias_strength = zeros(5, 1, 'double');
            % probability(step) = chance(step)/(2^32 - 1)
            p.bias_chance = zeros(5, p.time_steps, 'uint32'); 

            % Pre-allocate space for units
            p.unit_names = repmat({}, p.units, 1);            % Text name of unit
            p.unit_group_index = zeros(p.units, 1, 'uint16'); % Unit Index number within its group
            p.unit_threshold = zeros(p.units, 1, 'double');   % Each unit has its own threshold
            p.unit_bias_offset = zeros(p.units, 1, 'uint16'); % zero based bias id
            p.unit_pre_count = zeros(p.units, 1, 'uint16');   % Number of presynaptic units
            p.unit_post_count = zeros(p.units, 1, 'uint16');  % Number of postsynaptic units
            p.unit_pre_offset = zeros(p.units, 1, 'uint32');  % Zero based offset into weight_pre_sort
            p.unit_post_offset = zeros(p.units, 1, 'uint32'); % Zero based offset into wieght_post_sort
            p.unit_LFP_offset = zeros(p.units, 1, 'uint16');  % Zero based offset into LFP array
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

        end
        
        function initialize_stim_data(p, init_times, stim_replay_file)
            %INITIALIZE_STIM_DATA Pre-allocate for saving stim data
            %
            % Syntax:
            %   p.initialize_stim_data();
            %   p.initialize_stim_data(init_times, stim_replay_file);
            %
            % Inputs:
            %   init_times - (Default: [8.0, 8.7, 9.4] seconds) - This
            %                   gives the number of seconds until stimulus
            %                   test pulses begin happening for each
            %                   column. Values correspond to starts for
            %                   [Column A, Column B, Column C] respectively
            %   stim_replay_file - (Default: '') - This can be used to set
            %                   a "replay" stim ticks file that is useful
            %                   for "catch" trials using Spike-triggered or
            %                   Cycle-triggered stimulation with various
            %                   parameters.
            
            if nargin < 2
               init_times = [8.0, 8.7, 9.4]; 
            end
            
            if nargin < 3
               stim_replay_file = ''; 
            end
            
            % p.stim_source and p.target_times are used for paired pulse 
            %   stimulations and tetanic stimulation. 
            % p.stim_test_times are used for testing evoked potentials 
            %   with a column-wide stimulus during test sections.            
            p.stim_source_times = zeros(p.time_steps, 1, 'uint16');  % Default to no paired pulse or cycle triggered stims
            p.stim_target_times = zeros(p.time_steps, 1, 'uint16');
            p.stim_test_times =   zeros(p.time_steps, 1, 'uint16');
            p.stim_output_clockticks = zeros(0, 1, 'int32');         % Collect clocktick of each output stimulation.
            p.unit_spike_counts = zeros(1, p.units);                 % collect number of  output spikes for each unit
            p.stim_pulse_isi = floor(p.fs / p.stim_pulse_freq);      % stimulus pulse train interval.
            p.init_test_times = init_times;
            p.test_pulse_isi = floor(p.fs / p.test_pulse_freq); % Test pulse ISI in timesteps
            for istim = 0:p.test_pulse_train-1
                p.stim_test_times(floor(p.init_test_times(1) * p.fs + 1 + p.test_pulse_isi * istim)) = 1; % Col A test stimulus times start at 8 seconds
                p.stim_test_times(floor(p.init_test_times(2) * p.fs + 1 + p.test_pulse_isi * istim)) = 2; % Col B test stimulus times start at 8.7 seconds
                p.stim_test_times(floor(p.init_test_times(3) * p.fs + 1 + p.test_pulse_isi * istim)) = 3; % Col C test stimulus times start at 9.4 seconds
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
            p.stim_replay_file = stim_replay_file;
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
        end
        
    end
    
    methods (Hidden, Access=public)
        function pars = default_network_parameters(p)
           %DEFAULT_NETWORK_PARAMETERS Return default network parameters struct
           pars = struct(...
                'bias_correlated_fraction', p.bias_correlated_fraction, ...
                'bias_rate', p.bias_rate, ...
                'pair_uV', p.pair_uV_, ...
                'stim_excit_range', p.stim_excit_range_, ...
                'stim_inhib_range', p.stim_inhib_range_ ...
                );
        end
        
        function pars = default_time_factor_parameters(p)
           %DEFAULT_TIME_FACTOR_PARAMETERS Return default time factors struct
           pars = struct(...
                'PSP_slow_time_constant_ms', p.PSP_slow_time_constant_ms, ...
                'PSP_fast_time_constant_ms', p.PSP_fast_time_constant_ms, ...
                'STDP_strengthening_slow_time_constant_ms', p.STDP_strengthening_slow_time_constant_ms, ...
                'STDP_strengthening_fast_time_constant_ms', p.STDP_strengthening_fast_time_constant_ms, ...
                'STDP_weakening_slow_time_constant_ms', p.STDP_weakening_slow_time_constant_ms, ...
                'STDP_weakening_fast_time_constant_ms', p.STDP_weakening_fast_time_constant_ms, ...
                'training_rate', p.training_rate, ...
                'train_weakening', p.train_weakening, ...
                'train_pos_factor', p.train_pos_factor ...
            );
        end
        
        function pars = default_filter_parameters(p)
           %DEFAULT_FILTER_PARAMETERS Return default filtering parameters struct
           pars = struct(...
                'fs', p.fs, ...
                'axonal_delay_ms', p.axonal_delay_ms, ...
                'dendritic_delay_ms', p.dendritic_delay_ms, ...
                'output_delay_ms', p.output_delay_ms, ...
                'stim_delay_ms', p.stim_delay_ms, ...
                'stim_refractory_ms', p.stim_refractory_ms, ...
                'gamma_band', p.gamma_band, ...
                'beta_band', p.beta_band, ...
                'EMG_band', p.EMG_band ...
            );
        end
        
        function pars = default_normal_pdf_table_parameters(p)
            %DEFAULT_NORMAL_PDF_TABLE_PARAMETERS Return default parameters struct
            %
            % pars = p.default_normal_pdf_table_parameters();
           pars = struct( ...
               'correlated_bias_n_std', p.correlated_bias_n_std, ... % + / - 4 ==> 8
               'correlated_bias_std_ms', p.correlated_bias_std_ms ...
               );
        end
    end
    
    methods (Access=private)
        function init_network_(p)
            %INIT_NETWORK_ Private function to start network initialization
            %disp(['Estimated units: ' num2str(p.units) ', weights: ' num2str(p.weights)]);
            p.units = 0;
            p.weights = 0;
            p.biases = 1;
            % Caclulate PSP_factor for converting uV into weights for our
            % slow and fast exponential decay functions. This is based on the
            % peak value in the PSP curve calculated with the 0.1 ms timestep.
            % This keeps the areas under PSP curves with smaller timesteps more
            % comparable.
            slow_decay = 1 - 0.1 / p.PSP_slow_time_constant_ms; 
            fast_decay = 1 - 0.1 / p.PSP_fast_time_constant_ms;   
            timestep_index = p.PSP_comparison_window_ms; % Assumes PSP maximal value is within 100ms of start. 
            yslow = slow_decay .^ timestep_index;
            yfast = fast_decay .^ timestep_index;
            p.PSP_factor = 1 / max(yslow - yfast);  % Multiplier for PSP amplitudes to account for exponential decay function.
            p.max_PSP_value = p.max_strength * p.PSP_factor;
        end 
    end
    
    methods (Static)
        function s = update_fields_from_struct(s, pars)
            %UPDATE_FIELDS_FROM_STRUCT Update struct fields from second struct
            %
            % Syntax:
            %   s = SpikeNet.update_fields_from_struct(s, pars);
            %
            % Use this to update some fields or properties from an existing
            % properties (parameters) struct.
            
            f = fieldnames(pars);
            for iF = 1:numel(f)
                s.(f{iF}) = pars.(f{iF});
            end
        end
        
        function s = update_fields_from_args(s, varargin)
            %UPDATE_FIELDS_FROM_ARGS Update struct fields from input args
            %
            % Syntax:
            %   s = SpikeNet.update_fields_from_args(s, 'name', value, ...);
            %
            % Use this to update some fields or properties of an existing
            % struct or object using the varargin <'Name', value> syntax.
            %  Note:
            %   For a cell array in Matlab, using {:} puts all the
            %   cell elements in as subsequent arguments. So you can pass
            %   `varargin` of other methods to this method using the syntax
            %   s = SpikeNet.update_fields_from_args(s, varargin{:});
            %   When varargin is already an input argument to the calling
            %   function/method.
            if nargin < 2
                return;
            end
            if isstruct(varargin{1})
                s = varargin{1};
                varargin(1) = [];
            end
            for iV = 1:2:numel(varargin)
                s.(varargin{iV}) = varargin{iV+1};
            end
        end
    end
end

