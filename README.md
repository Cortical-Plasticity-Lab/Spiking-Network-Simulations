# Spiking Network Simulations #

_Sections_
* **[Code Summary](#simulation-scripts)**
* **[Running Simulations](#use)**
* **[Troubleshooting](#setup)**

## Developer Notes ##
* 2021-04-14: Forked repository.
* *This code has been adapted from the original source code.*
  *  **Original code** is provided by **[Larry Shupe](https://github.com/lshupe/spikenet)**. 
  * **Associated original manuscript** has been published on **[eNeuro](https://www.eneuro.org/content/8/2/ENEURO.0333-20.2021)**.
* *This README was adapted from the original `spikenet_readme.txt`, with some modifications.*

## Contents ##

### Simulation Scripts ###

These MATLAB `script` files that run simulations with different pre-configured parameters.

* `spikenet50.m`  -- Run the base simulation without conditioning.
* `spikenet_spike_trig.m`  -- Configured for spike triggered conditioning
* `spikenet_paired_pulse.m`  -- Configured for paired pulse conditioning
* `spikenet_cycle_trig.m`  -- Configured for cycle triggered conditioning
* `spikenet_tetanic.m`  -- Configured for tetanic conditioning
* `spikenet_emg_trig`  -- Configured for EMG triggered conditioning
* `spikenet_gamma_trig`  -- Configured for gamma filtered LFP triggered conditioning
* `spikenet_FICspike_trig.m`  -- Configured for spike triggered conditioning with Fixed Intracolumn Connections (FIC).

### Helpers ###

* `spikenet50mex.c`  -- The `mex` function that runs the network in blocks of 10 seconds
* `spikenet50mex.mexw64`  -- Precompiled MATLAB `mex` function, which can be used to re-compile `spikenet50mex.c` in case you are running a different OS or version of MATLAB.
* `spikenet_wfig.m`  -- Creates a weight figure from a simulation’s saved `param.mat` file.

## Use ##

Once the simulation files have been installed, open MATLAB and change the current 
directory to the folder containing these files (or add the files to the default MATLAB `path` as described **[above](#setting-the-matlab-path)**).  

### Running a simulation without conditioning ###

1. Run `spikenet_no_cond.m`.
   * This script will take up to an hour to run, but will update a graph 
     for every 10 seconds of simulated time.  This graph is used to monitor the progress of the simulation.  It displays average activity of excitatory and output units for the 10 second time blocks for the current training period, averages of filtered local field potentials (LFP), and the current distribution of trainable connection strengths.  
   * These averages reset at the beginning of each of the four training 
     periods: 
     * No Conditioning, 
     * Preconditioning testing, 
     * Conditioning, and 
     * Post-Conditioning testing 
   * Connections between units are modified by the Spike-Timing Dependent 
     Plasticity (STDP) rule during the Preconditioning and Conditioning periods, but 
     connection strengths are held static during the testing periods to build histograms and conduct tests.  

### Data structure and location ###

#### Data Location ####

> `C:/data/directory`

#### Figure Name Schema ####

All files ending in `.fig` are MATLAB `figure` files that can be opened directly in the MATLAB base workspace to reproduce the vectorized figures (and associated data).

*Note*: in each of the following filenames, (`<n>` is the the simulation trial number for each 10 second block).

* `averages_cond_<n>.fig` – The  monitoring averages saved during the conditioning period.
* `averages_nocond_<n>.fig` – The monitoring averages saved during the no-conditioning period. 
* `averages_pre_<n>.fig` – The monitoring averages saved during the pre-conditioning period.
* `averages_post_<n>.fig` – The monitoring averages saved during the post-conditioning period.
* `Cond_Trig_PSTH-<n>.fig` – Spike triggered averages and histograms aligned with conditioning triggers.
* `emg_spike_response.fig` – Stimulus triggered averages of the simulated EMG.
* `lfp_spike_response.fig` – Spike-triggered LFP averages.
* `lfp_stim_response.fig` – Stimulus triggered LFP averages (used for A->B evoked potential changes).
* `mean_Evoked_Potentials.fig` – Tracks MPI from the lfp_stim_response.fig when multiple repetitions through the training sequence are used (when p.ntrials_repeat > 1)
  output_spikes.mat – Contains spike times (in 0.1 ms time steps) for unit firings during the simulation.
* `param.mat` – Contains a copy of the simulations parameters and resulting connection strengths.
* `progress.txt` – A copy of the console output during the simulation.
* `PSTH.fig` – Shows correlated activity between unit populations
* `Raster-<n>.fig` -- Contains a dot raster for 2 seconds of activity and power spectrums for LFP.
* `weight_distribution.fig` – Distribution of weights between unit columns.  
  * This is useful for checking the efficacy of any arbitrary conditioning method of interest in greater detail.
* `weight_matrix.fig` – Displays final connection strengths.
  * Use `spikenet_wfig.m` on the `param.mat` file to plot a formatted graph with labels.

#### Other Saved Outputs ####

These files are primarily used for validation and checking consistency/reproducibility of the simulations that were conducted.

* `spiknet.m` – Contains a copy of the `function` that runs this simulation.
  * This is useful to validate/verify that the correct parameters were implemented for a given simulation iteration. 
* `spikenet50mex.c` – contains a copy of the `mex` function.
* `spikenet50mex.mexw64` – Contains a copy of the compiled `mex` function.

## Experiments ##

### Conditioning Methods ###

*Contents*

* **[No Conditioning](#no-conditioning)**
* **[Spike-triggered Stimulation](#spike-triggered-stimulation)**
* **[Paired-Pulse Stimulation](#paired-pulse-stimulation)**
* **[Cycle-triggered Stimulation](#cycle-triggered-stimulation)**
* **[Tetanic Conditioning](#tetanic-conditioning)**
* **[EMG-triggered Stimulation](#emg-triggered-stimulation)**
* **[Gamma-triggered Stimulation](#gamma-triggered-stimulation)**
* **[Spike-triggered Stimulation (Fixed Inter-columnar Connections; FIC)](#spike-triggered-fic-stimulation)**

#### No Conditioning ####

The base (non-conditioning) simulation provides information on network progression in the absence of any conditioning stimulus intervention.  This can be used to see how much random variation there is by running the network for different periods of time or by using different initial conditions.  The base simulation file `spikenet_no_cond.m` contains a standard set of parameters. **These parameters are used by the other simulations.** 

To trouble-shoot parameter configurations, start here!

[*back to list of stimulation methods*](#conditioning-methods)

#### Spike-triggered Stimulation #### 

This conditioning method delivers a stimulus to units in **Column B** each time that the first excitatory unit in **Column A** (`Ae1`) fires. After a given stimulus, any subsequent stimuli are "blanked" (do not occur) for the duration of the "blanking period" parameter. This prevents stimulus-locked feedback loops from occurring, which could be deleterious *in vivo*. Since the timing of the stimuli depend on detected activity in the network, this is referred to as a “closed loop” stimulation method.

##### Notable parameters #####
```matlab
% % % Key to running this simulation! % % %
p.conditioning_method = 1;  	% Selects spike triggered stimulation

% % Also important % %
p.stim_delay_ms = 10;    	% Millisecond delay between spike on Ae1 and
				% 	stimulus on Column B.
p.stim_pulse_train = 1;  	% Can be 2 or 3 for stimulus trains of 2 or 3
				% 	pulses per train.
p.stim_uV = 2000;        	% Size of conditioning stimulus 
   				%	(0 can be used as a sham stimulus)
p.stim_refractory_ms = 10;  	% Refractory period on delivered stimulation.
```

[*back to list of stimulation methods*](#conditioning-methods)

#### Paired-Pulse Stimulation #### 

This conditioning method delivers a stimulus to units in **Column A** followed by a delayed stimulus that is delivered to units of **Column B**.  These stimuli are not tied to any activity in the network. This is therefore an “open loop” stimulation method (with respect to the "natural" activity of the network).

##### Notable parameters #####

```matlab
% % % Key to running this simulation! % % %
p.conditioning_method = 2;  		% Selects paired pulse stimulation

% % Also important % %
p.conditioning_secs = 7; 		% 1..10 to change number of paired pulses in each 10 second time block.
p.stim_delay_ms = 10;    		% Delay between the paired pulses
					%	Units: milliseconds
p.stim_pulse_train = 1;  		% Number of stimulus pulses per train. 
					%	Although the blanking period "debounces" additional stimulus trains
                                	%	increasing this count greater than one will still deliver multiple
                                	%	stimuli each time a stimulus is "triggered."
p.stim_uV = 2000;        		% Size of conditioning stimulus (0 can be used as a sham stimulus)
					%	Units: micro-volts
p.pair_uV = 2000;        		% Size of the paired pulse simulation (the second pulse)
					%	Units: micro-volts
```

#### Cycle-triggered Stimulation ####

A "closed loop" conditioning method that detects oscillations in the spike rate of **Column B** in order to trigger a stimulus train targeting units in **Column A**.

##### Notable parameters #####

```matlab
% % % Key to running this simulation! % % %
p.conditioning_method = 3;  	% Selects cycle triggered stimulation

% % Also important % %
p.bias_modulation_rate_B = 20; 	% Cycles per second modulation rate on Column B. 0 = no modulation. 
				% 	Units: Hz
p.bias_modulation_amount = .2; 	% Size of modulation. 0.2 = +/- 20% of normal rate.
p.stim_pulse_train = 1;        	% Number of stimulus pulses per train. 
				%	Although the blanking period "debounces" additional stimulus trains
                                %	increasing this count greater than one will still deliver multiple
                                %	stimuli each time a stimulus is "triggered."
p.stim_uV = 2000;              	% Size of conditioning stimulus (0 can be used as a sham stimulus)
				%	Units: micro-volts
```

[*back to list of stimulation methods*](#conditioning-methods)

#### Tetanic Conditioning ####

An open loop conditioning method that applies random stimulation on the units of a specific column.  This is an important control for other stimulation methods, because it is completely "open loop" with respect to both network activity and any pairing of **Column A** with **Column B** (unlike occurs in **[Paired-Pulse Stimulation](#paired-pulse-stimulation)**)

##### Notable parameters #####

```matlab
% % % Key to running this simulation! % % %
p.conditioning_method = 5; 	% Tetanic stimulation: 
				%   4 => Stimulation targets Column A units
				%   5 => Stimulation targets Column B units
				%   6 => Stimulation targets Column C units

% % Also important % %
p.tetanic_freq = 10;       	% Stimulation intensity (rate parameter).
				% 	This sets the intensity of the exponentially distributed 
				%	inter-stimulus intervals in the tetanic conditioning regime.
				%		Units: stimuli/sec
p.stim_refractory_ms = 10; 	% Refractory period.
				%	Sets the minimum number of milliseconds between tetanic stimuli. 
				%	This effectively reduces the average rate of stimulation, 
				%	since it prevents realization of a true exponential distribution. 
				%	It is included to reproduce the tetanic stimulation refractory period
                                %	from in vivo protocols, which is included as a 
                                %	precautionary safety measure.
				%		Units: milliseconds
p.stim_pulse_train = 1;    	% Number of stimulus pulses per train. 
				%	Although the blanking period "debounces" additional stimulus trains
                                %	increasing this count greater than one will still deliver multiple
                                %	stimuli each time a stimulus is "triggered."
p.stim_uV = 2000;          	% Size of conditioning stimulus (0 can be used as a sham stimulus)
				%	Units: micro-volts
p.bias_modulation_amount = 0; 	% Normally 0.  Can put a 0.2 to 0.4 step bias spike increase on Column B.
				%	Units: arbitrary (probability)
```

[*back to list of stimulation methods*](#conditioning-methods)

#### EMG-triggered Stimulation ####

A closed loop conditioning method that applies stimuli to units in **Column B** when the output of the transfer function relating activity of units in **Column A** to the generation of electromyographic signals in some target muscle(s) crosses a specified threshold.

##### Notable parameters #####

```matlab
p.conditioning_method = 7;  	% Select EMG triggered stimulation
p.lfp_detect_level = 1000;  	% EMG detection level
				%	Units: arbitrary (see transfer function for EMG output)
p.stim_delay_ms = 0;        	% Delay between EMG level detection and stimulation of Column B. 
				%	(Note: The transfer function generating synthetic EMG signals 
				%				intrinsically delays the output EMG by 10-ms due to the 
				%				default synaptic delay parameters between Ae and Ao units)
				%	Units: milliseconds
p.stim_pulse_train = 1;  	% Number of stimulus pulses per train. 
                            	%	Although the blanking period "debounces" additional stimulus trains
                            	%	increasing this count greater than one will still deliver multiple
                            	%	stimuli each time a stimulus is "triggered."
p.stim_uV = 2000;        	% Size of conditioning stimulus (0 can be used as a sham stimulus)
				%	Units: micro-volts
```

[*back to list of stimulation methods*](#conditioning-methods)

#### Gamma-triggered Stimulation ####

A "closed loop" conditioning method that applies stimuli to units in **Column B** when the band-passed local field potential (LFP; generated via a transfer function relating the activity of units in **Column A** to a field potential source capable of generating a time-varying electrical field) crosses an arbitrary threshold.

##### Notable parameters #####

```matlab
% % % Key to running this simulation! % % %
p.conditioning_method = 8;  	% Select gamma triggered stimulation

% % Also important % %
p.stim_delay_ms = 0;          	% Millisecond delay between gamma threshold-level detection and 
				%	stimulation delivery to Column B. 
				%		(use value 8+ for rising level detection)
				%			Units: milliseconds
p.lfp_detect_level = -20000;  	% Filtered LFPA detection level (20000 for rising level, -20000 for falling)
p.gamma_band = [50 80];       	% Bandwidth that selects for `gamma` from the broadband local field
				%	potential (LFP) signal. (As an example, [50 80] specifies cutoff
				%	frequencies of 50Hz and 80Hz for the bandpass filter).
				%			Units: Hz
p.stim_pulse_train = 1;  	% Number of stimulus pulses per train. 
                            	%	Although the blanking period "debounces" additional stimulus trains
                            	%	increasing this count greater than one will still deliver multiple
                            	%	stimuli each time a stimulus is "triggered."
p.stim_uV = 2000;        	% Size of conditioning stimulus (0 can be used as a sham stimulus)
				%	Units: micro-volts
```

[*back to list of stimulation methods*](#conditioning-methods)

#### Spike-triggered FIC Stimulation #### 

Similar to normal spike-triggered conditioning ("closed loop"), except that fixed intercolumnar connections (FIC) are used in order to simulate the role of the correlated bias inputs to provide correlated activity within each column.

```matlab
% % % Key to running this simulation! % % %
p.conditioning_type = 1; 			% Note: this is the same as the default "no conditioning" condition, but
						%	uses altered parameters.

% % Also important % %
p.bias_corrrelated_fraction = 0.0
p.initmin = 200;  					% Lower bound on initial PSP amplitude for any connections
p.initmax = 300; 					% Upper bound on initial PSP amplitude for any connections
p.epsp2excit_incol = [p.initmin p.initmax 1 0];     	% EPSP for connections to excitatory units within a column
p.epsp2inhib_incol = [p.initmin p.initmax 1 0];     	% EPSP for connections to inhibitory units within a column
p.ipsp2excit_incol = [-p.initmax -p.initmin 1 0];   	% IPSP for connections to excitatory units within a column
p.ipsp2inhib_incol = [-p.initmax -p.initmin 1 0];   	% IPSP for connections to inhibitory units within a column
p.stim_delay_ms = 10;    				% Delay between spike on Ae1 and stimulus on Column B.
							%	Units: milliseconds
p.stim_pulse_train = 1;  	% Number of stimulus pulses per train. 
                        	%	Although the blanking period "debounces" additional stimulus trains
				%	increasing this count greater than one will still deliver multiple
                            	%	stimuli each time a stimulus is "triggered."
p.stim_uV = 2000;        	% Size of conditioning stimulus (0 can be used as a sham stimulus)
				%	Units: micro-volts
p.stim_refractory_ms = 10;  	% Refractory period on delivered stimulation.
				%	Units: milliseconds
```

[*back to list of stimulation methods*](#conditioning-methods)

### General Comments ###

Most parameters are set **within the first 200 lines** of the **[main simulation script](#no-conditioning)**.  

#### Other Parameters ####

The following is a list of parameters that have not yet been described. Since they have been modified in previous simulations, it is likely that they will be altered again in future simulations.

The following parameters are used to control the length of each training phase. Larger values for `p.ntrials_off` will allow the network to "settle" for longer during the pre-conditioning phase. 

**Note about initial phase (conditioning off) duration.**

Larger values of `p.ntrials_off` allow better convergence to "natural" synaptic connection values given the network topology. However, convergence comes at the expense of longer simulation times!

```matlab
p.ntrials_off = 50;     % Number of initial trials with conditioning off (STDP training on).
p.ntrials_on = 50;      % Number of following trials with conditioning on (STDP training on).
p.ntrials_test = 50;    % Number of trials for summary figures (both conditioning and STDP training off)
p.ntrials_repeat = 1;   % Number of repetitions through the off/test/on/test cycle.
```

The firing threshold for units can be increased to lower overall firing rates, but must be balanced against the external bias input rates.  The correlated external biases are only correlated within each individual column.  

Increasing the correlated fraction will tend to increase firing rates and also improve conditioning effects that depend on intercolumn correlated activity such as spike-triggered stimulation. 

Output units only receive uncorrelated bias spikes, the rate of which can be adjusted independently.

```matlab
p.threshold = 5000;   				% Threshold for all cortical units.
						%	Units: micro-volts
p.bias_value = 350;   				% Strenth of bias potientials.
						%	Units: micro-volts
p.bias_rate = 1800;   				% Bias of each column unit. 
						%	Consider smaller values if biases are modulated.
						%	Units: spikes/second
p.bias_correlated_fraction = .30; 		% Fraction of correlated bias spikes
p.correlated_bias_std_ms = 3;      		% Standard deviation of normally distributed (correlated) spiking bias.
						%	Units: milliseconds
p.out_bias_rate = 2000;            		% Uncorrelated bias spikes for output units.
						%	Units: spikes/second
```

The training rate and weakening factors can be altered to change the how fast connection strengths change.  Decreasing the training_rate usually means the network will need to train longer, but if the value is too high then there will be more variability in weights and it may seem like the network doesn’t settle down.  Increasing the weakening factor controls the effectiveness of the 
weakening side (negative post-presynaptic spike time differences) of the STDP curve relative to the strengthening side (positive post-presynaptic spike time differences).  

Decreasing the weakening factor can increase the overall activity of the network, but at a certain point all weights will tend to increase to the maximum allowed value.  Typically we want to have the weakening factor large enough to prevent that but not so large that connections in the network cannot be conditioned.  This will also be affected by the shape of the curve which is 
controlled by fast and slow decay constants for both halves of the STDP function.

```matlab
p.training_rate = 100;    		% Train rate factor for spike-timing-dependent plasticity. 
					% 	This factor accomodates both the strengthing (positive) and 
					%	weakening (negative) sides of the SDTP rule, 
					%	so it is always specified as a positive value.
p.train_weakening = 0.55; 		% Relative amplitude for the weakening side of the SDTP curve.
```

The network topology can be changed by altering the column connection definitions.

Initial weights are from a uniform distribution between a minimum and maximum value.  

The connection probability for excitatory units is usually 1/6, but can be increased to provide increased network activity.  Inhibitory connections are only within column and have a larger connection probability to help balance the number of inhibitory vs excitatory connections in the network.  

The random seed allows the use of the same set of pseudo-random connections for many different simulations.  The initial network configuration *can have a significant effect on the outcome*, **so being able to reuse a network configuration is important**. In practice, simulations should be run **multiple times with *different* random seeds**, then average together the results from those random seeds. However, for reproducibility, it is important to keep track of the random seeds that were used!

```matlab
% Initial post-spike potential (PSP) ranges in uV.  
p.initmin = 100;		% Initial minimum PSP value (micro-volts)
p.initmax = 300; 		% Initial maximum PSP value (micro-volts)

% All array parameters follow the format:
%	[microvolts_min microvolts_max probability_of_connection enable_training]

% EPSPs for connections to excitatory units within a column
p.epsp2excit_incol = [p.initmin p.initmax 1/6 1];     

% EPSPs for connections to inhibitory units within a column
p.epsp2inhib_incol = [p.initmin p.initmax 1/6 1];     

 % EPSPs for connections to excitatory units in adjacent columns
p.epsp2excit_outcol = [p.initmin p.initmax 1/6 1];  

% EPSPs for connections to inhibitory units in adjacent columns
p.epsp2inhib_outcol = [p.initmin p.initmax 1/6 1];    

% IPSPs for connections to excitatory units within a column
p.ipsp2excit_incol = [-p.initmax -p.initmin 1/3 1];   

% IPSPs for connections to inhibitory units within a column
p.ipsp2inhib_incol = [-p.initmax -p.initmin 1/3 1];   

% EPSPs for connections to output units.
p.epsp2output_incol = [350 350 1/3 0];      
% Note: These values are **not** usually trained. 
%	Negative connection chance is interpreted by the simulations as:
%	 > "Connect to the first p * N units."
%	Therefore, it is a way to set a FIXED set of connections, rather than chance connections.

p.random_seed = 1; 	% Fixed random number seed. 
			%	This allows us to use the same network 
			%	topology (connections) for different conditioning methods.
			%	
			%	The random seed should be altered at different stages of the simulations, 
			%	depending on which factors are important to test and where the variability arises.
```

## Setup

### Software Requirements ###

* Windows 10 (64-bit)

  > The source code for the integrate-and-fire (IF) neural network model is written in Matlab and C for the 64-bit version running under Windows 10. 

* Matlab (R2017a+)

  > The C code is written as a Matlab `mex` function, which has been pre-compiled for current (R2017a+) 64-bit versions of Matlab for PC. 

 ### Troubleshooting ###

#### Parameter configuration ####

The base simulation file `spikenet_no_cond.m` contains a standard set of parameters. **These parameters are used by the other simulations.** To trouble-shoot parameter configurations, start here!

#### Re-compiling `mex` function ####

For older versions of Matlab, or for other operating systems such as Mac OSX, it is necessary to re-compile the `mex` file, which allows the code to run efficiently.

> This can be accomplished by using the `mex` command from within the MATLAB Command Window interface, specifying the correct C compiler. 

*Note:* the system must have a valid compiler configured, which is then set up in MATLAB separately from the `mex` command. For more information, please see the MATLAB `mex` documentation, located **[here](https://www.mathworks.com/help/matlab/ref/mex.html)**.

#### Setting the MATLAB `path` ####

> To install the source code, copy all files to a folder on the MATLAB `path`.
>
> * This is usually the MATLAB folder within the user’s Documents folder (e.g. `C:/<user>/Documents/Matlab/` on a PC).

#### About MATLAB `script` and `function` files ####

(**For beginners**): MATLAB `.m` files that are not part of an object-oriented structure are either `script` files or `function` files. 

* Unlike MATLAB a `function`, a MATLAB `script` does not require any input arguments and therefore can be treated like an executable program in MATLAB. 

* To run a MATLAB `script`, type the name (without the file extension) of the `script` into the `Command Window` directly in the MATLAB editor, and press enter.

* All of the `.m` files in **this project** are MATLAB `function` files. A MATLAB `function` starts with the syntax `function <filename>(<arguments>)`.

  ```matlab
  % Example: if the file is `spikenet50.m` 
  % function spikenet50()
  ```

  A MATLAB `function` that does not take any arguments can be called in the same way as a MATLAB `script`, and the only difference is that all of the `function` variables will only exist in the `function` workspace, as opposed to the base MATLAB workspace as they would if the `function` were instead a `script`. These functions can be invoked in the same way as a MATLAB `script`:

  ```matlab
  >> spikenet_no_cond;
  ```

  Alternatively, the `script` can be run by opening it in the MATLAB editor and clicking the green `Run` button at the top of the `Editor` tab. 

  *Note:* using MATLAB `function` files instead of `script` files makes it less-likely that you will use incorrect parameters from one simulation while running a different simulation and is therefore better practice, although the use of `script` files can make debugging and development go faster since the `script` is effectively always running in "debug" mode with variables appearing in the Editor's Base workspace.
