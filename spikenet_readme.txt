Source Code for the Integrate and Fire Spiking Neural Network
Available at https://github.com/lshupe/spikenet

The source code for the IF neural network mode is written in Matlab and C for the
64-bit version running under Windows 10.  The C code is written as a Matlab mex 
function which has been precompiled for current (2017 to 2020) 64-bit versions of 
Matlab for PC.  This can be compiled with the mex command from within Matlab with 
the C compiler add-in.

To install the source code, copy all files to a folder on the Matlab Path.
This is usually the Matlab folder within the user’s Documents folder.  
For example C:/<user>/Documents/Matlab/ on a PC.

Contents:
spikenet50.m  -- This is the base simulation file with conditioning turned off.
spikenet_spike_trig.m  -- Configured for spike triggered conditioning
spikenet_paired_pulse.m  -- Configured for paired pulse conditioning
spikenet_cycle_trig.m  -- Configured for cycle triggered conditioning
spikenet_tetanic.m  -- Configured for tetanic conditioning
spikenet_emg_trig  -- Configured for EMG triggered conditioning
spikenet_gamma_trig  -- Configured for gamma filtered LFP triggered conditioning
spikenet50mex.c  -- The mex function that runs the network in blocks of 10 seconds
spikenet50mex.mexw64  -- Precompiled Matlab mex function
spikenet_wfig.m  -- Creates a weight figure from a simulation’s saved param.mat file.
spikenet_readme.txt  -- Information about running the network simulations.

Running a simulation
Once the simulation files have been installed, open Matlab and change the current 
directory to the folder containing these files.  Run the spikenet_no_cond.m script 
directly from the command line or open the script and run it from the editor tab 
run button.  This script will take up to an hour to run, but will update a graph 
for every 10 seconds of simulated time.  This graph is used to monitor the progress 
of the simulation.  It displays average activity of excitatory and output units for 
the 10 second time blocks for the current training period, averages of filtered 
local field potentials (LFP), and the current distribution of trainable connection 
strengths.   These averages reset at the beginning of each of the four training 
periods: No Conditioning, Preconditioning testing, Conditioning, and Post-Conditioning
testing.  Connections between units are modified by the Spike-Timing Dependent 
Plasticity (STDP) rule during the Preconditioning and Conditioning periods, but 
connection strengths are held static during the testing periods to build histograms 
and make tests.  

Simulation results are stored in a folder in the c:/data/ directory
(<n> is the the simulation trial number for each 10 second block):
averages_cond_<n>.fig – The monitoring averages saved during the conditioning period
averages_nocond_<n>.fig – The monitoring averages saved during the no-conditioning period
averages_pre_<n>.fig – The monitoring averages saved during the pre-conditioning period
averages_post_<n>.fig – The monitoring averages saved during the post-conditioning period
Cond_Trig_PSTH-<n>.fig – Spike triggered averages and histograms aligned with conditioning triggers
emg_spike_response.fig – Stimulus triggered averages of the simulated EMG
lfp_spike_response.fig – Snit triggered LFP averages
lfp_stim_response.fig – Stimulus triggered LFP averages (used for A->B evoked potential changes)
mean_Evoked_Potentials.fig – Tracks MPI from the lfp_stim_response.fig when multiple repetitions through the training sequence are used (when p.ntrials_repeat > 1)
output_spikes.mat – Contains spike times (in 0.1 ms time steps) for unit firings during the simulation
param.mat  – Contains a copy of the simulations parameters and resulting connection strengths
progress.txt – A copy of the console output during the simulation.
PSTH.fig – Shows correlated activity between unit populations
Raster-<n> -- Contains a dot raster for 2 seconds of activity and power spectrums for LFP.
spiknet.m – Contains a copy of the script that runs this simulation
spikenet50mex.c – contains a copy of the mex function
spikenet50mex.mexw64 – Contains a copy of the compiled mex function
weight_distribution.fig – Distribution of weights between unit columns.  Used for checking efficacy of the conditioning method
weight_matrix.fig – Displays final connection strengths.  Use spikenet_wfig.m on the param.mat file to plot a formatted graph with labels.


Conditioning Methods
The base no-conditioning simulation provides information on network progression 
in absence of any conditioning stimulus intervention.  This can be used to see 
how much random variation there is by running the network for different periods 
of time or by using different initial conditions.  The base simulation file 
spikenet_no_cond.m contains a standard set of parameters that the other simulation 
files modify.

Spike triggered stimulation:  This conditioning method delivers a stimulus to units 
in Column B when ever the first excitatory unit in Column A (Ae1) fires.  Since the 
timing of the stimuli depend on detected activity in the network, this is referred 
to as a “closed loop” stimulation method.
Notable parameters:
   p.conditioning_method = 1;  % Select spike triggered stimulation
   p.stim_delay_ms = 10;    % Millisecond delay between spike on Ae1 and stimulus on Column B.
   p.stim_pulse_train = 1;  % Can be 2 or 3 for stimulus trains of 2 or 3 pulses per train.
   p.stim_uV = 2000;        % Size of conditioning stimulus (0 can be used as a sham stimulus)
   p.stim_refractory_ms = 10;  % Refractory period on delivered stimulation.


Paired pulse stimulation:  This conditioning method delivers a stimulus to units 
in Column A followed by a delayed stimulus on Column B.  Since these stimuli are 
not relative to any activity in the network, this is an “open loop” stimulation method.
Notable parameters:
   p.conditioning_method = 2;  % Selects paired pulse stimulation
   p.conditioning_secs = 7; % 1..10 to change number of paired pulses in each 10 second time block.
   p.stim_delay_ms = 10;    % Millisecond delay between the paired pulses
   p.stim_pulse_train = 1;  % Can be 2 or 3 for paired pulse trains of 2 or 3 pulses per train.
   p.stim_uV = 2000;        % Size of conditioning stimulus (0 can be used as a sham stimulus)
   p.pair_uV = 2000;        % Size of the paired pulse simulation (the second pulse)


Cycle triggered stimulation:  A closed loop conditioning method that detects 
oscillations imposed on Column B to trigger a stimulation on Column A.
Notable parameters:
   p.conditioning_method = 3;  % Selects cycle triggered stimulation
   p.bias_modulation_rate_B = 20; % Cycles per second modulation rate on Column B. 0 = no modulation. 20 for 20Hz sine wave
   p.bias_modulation_amount = .2; % Size of modulation. 0.2 = +/- 20% of normal rate.
   p.stim_pulse_train = 1;        % Can be 2 or 3 for stimulus trains of 2 or 3 pulses per train.
   p.stim_uV = 2000;              % Size of conditioning stimulus (0 can be used as a sham stimulus)


Tetanic conditioning:  An open loop conditioning method that applies random 
stimulation on the units of a specific column.  This is an important control 
for other stimulation methods.
Notable parameters:
   p.conditioning_method = 5; % Tetanic stimulation: 4 = Column A, 5 = Column B, 6 = Colunm C
   p.tetanic_freq = 10;       % Max stims/sec for exponentially distributed tetanic conditioning 
   p.stim_refractory_ms = 10; % Refractory period will reduce the actual stimulation frequency.
   p.stim_pulse_train = 1;    % Can be 2 or 3 for stimulus trains of 2 or 3 pulses per train.
   p.stim_uV = 2000;          % Size of conditioning stimulus (0 can be used as a sham stimulus)
   p.bias_modulation_amount = 0; % Normally 0.  Can put a 0.2 to 0.4 step bias spike increase on Column B.


EMG triggered conditioning:  A closed loop conditioning method that applies 
stimuli on Column B when the synthetic EMG level of Column A crosses a threshold.
Notable parameters:
   p.conditioning_method = 7;  % Select EMG triggered stimulation
   p.lfp_detect_level = 1000;  % EMG detection level
   p.stim_delay_ms = 0;        % Millisecond delay between EMG level
       detection and stimulus on Column B. (EMGs are already delayed 10 ms
       by synaptic delay between Ae and Ao units)
   p.stim_pulse_train = 1;  % Can be 2 or 3 for stimulus trains of 2 or 3 pulses per train.
   p.stim_uV = 2000;        % Size of conditioning stimulus (0 can be used as a sham stimulus)


Gamma triggered conditioning:  A closed loop conditioning method that applies 
stimuli to Column B when the band-passed LFP of Column A crosses a threshold.
Notable parameters:
   p.conditioning_method = 8;  % Select gamma triggered stimulation
   p.stim_delay_ms = 0;          % Millisecond delay between gamma level
       detection and stimulus on Column B. (use value 8+ for rising level detection)
   p.lfp_detect_level = -20000;  % Filtered LFPA detection level (20000 for rising level, -20000 for falling)
   p.gamma_band = [50 80];       % Bandwidth for conditioning_type 8 in Hz (for example [50 80] for 50Hz to 80Hz bandpass filter).
   p.stim_pulse_train = 1;  % Can be 2 or 3 for stimulus trains of 2 or 3 pulses per train.
   p.stim_uV = 2000;        % Size of conditioning stimulus (0 can be used as a sham stimulus)


Most parameters are set within the first 200 lines of the main simulation script.  
The following is a list of parameters that have not yet been described but have 
been modified in some previously run simulations and so are likely to be altered 
in future simulations.

The following parameters are used to control the length of each training phase. 
Larger values for p.ntrials_off will allow the network to settle for longer during 
the preconditioning phase at the expense of longer simulation times.

p.ntrials_off = 50;     % Number of initial trials with conditioning off (STDP training on).
p.ntrials_on = 50;      % Number of following trials with conditioning on (STDP training on).
p.ntrials_test = 50;    % Number of trials for summary figures (both conditioning and STDP training off)
p.ntrials_repeat = 1;   % Number of repetitions through the off/test/on/test cycle.


The firing threshold for units can be increased to lower overall firing rates, 
but must be balanced against the external bias input rates.  The correlated 
external biases are only correlated within each individual column.  Increasing 
the correlated fraction will tend to increase firing rates and also improve 
conditioning effects that depend on intercolumn correlated activity such as 
spike-triggered stimulation.  Output units only receive uncorrelated bias spikes, 
the rate of which can be adjusted independently.

p.threshold = 5000;   % uV threshold for all cortical units
p.bias_value = 350;   % uV strenth of bias potientials
p.bias_rate = 1800;   % Bias spikes per second to each column unit. Consider smaller values if biases are modulated.
p.bias_corrrelated_fraction = .30; % Fraction of correlated bias spikes
p.correlated_bias_std_ms = 3;      % Standard deviation of the normally distributed (correlated) bias spikes (in milliseconds)
p.out_bias_rate = 2000;            % Number of uncorrelated bias spikes per second for output units.


The training rate and weakening factors can be altered to change the how fast 
connection strengths change.  Decreasing the training_rate usually means the 
network will need to train longer, but if the value is too high then there will 
be more variability in weights and it may seem like the network doesn’t settle 
down.  Increasing the weakening factor controls the effectiveness of the 
weakening side (negative post-presynaptic spike time differences) of the STDP 
curve relative to the strengthening side (positive post-presynaptic spike time 
differences).  Decreasing the weakening factor can increase the overall activity 
of the network, but at a certain point all weights will tend to increase to the 
maximum allowed value.  Typically we want to have the weaking factor large enough 
to prevent that but not so large that connections in the network cannot be 
conditioned.  This will also be affected by the shape of the curve which is 
controlled by fast and slow decay constants for both halves of the STDP function.

p.training_rate = 100;    % Train rate factor for both strengthing (pos) and weakening (neg) sides of the SDTP rule.
p.train_weakening = 0.55; % Relative amplitude for the weakening side of the SDTP curve.


The network topology can be changed by altering the column connection definitions.  
Initial weights are from a uniform distribution between a minimum and maximum value.  
The connection probability for excitatory units is usually 1/6, but can be increased 
to provide increased network activity.  Inhibitory connections are only within column 
and have a larger connection probability to help balance the number of inhibitory vs 
excitatory connections in the network.  The random seed allows the use of the same 
set of pseudo-random connections for many different simulations.  The initial network 
configuration can have a significant effect on the outcome, so being able to reuse a 
network configuration is important.

% Initial PSP ranges in uV.  Arrays are [uVmin uVmax pChanceForConnection enableTraining]
p.initmin = 100;
p.initmax = 300;
p.epsp2excit_incol = [p.initmin p.initmax 1/6 1];     % Epsp for connections to excitatory units within a column
p.epsp2inhib_incol = [p.initmin p.initmax 1/6 1];     % Epsp for connections to inhibitory units within a column
p.epsp2excit_outcol = [p.initmin p.initmax 1/6 1];    % Epsp for connections to excitatory units in adjacent columns
p.epsp2inhib_outcol = [p.initmin p.initmax 1/6 1];    % Epsp for connections to inhibitory units in adjacent columns
p.ipsp2excit_incol = [-p.initmax -p.initmin 1/3 1];   % Ipsp for connections to excitatory units within a column
p.ipsp2inhib_incol = [-p.initmax -p.initmin 1/3 1];   % Ipsp for connections to inhibitory units within a column
p.epsp2output_incol = [350 350 1/3 0];    % Epsp for connections to output units.  These are not usually trained. Negative connection chance means first p * N units rather than random chance.
p.random_seed = 1; % Fixed random number seed so we can use the same network topology for different conditioning methods.


