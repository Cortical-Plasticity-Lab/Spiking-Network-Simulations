/*******************************************************************************
% spikenet50mex.c	
% Larry Shupe
% Aug 14, 2017
% Matlab .mex file 
% Copyright (C) 2015 University of Washington Regional Primate Research Center.
% All rights reserved.
%
% The Matlab calling syntax is:
%   [uindex, tstamp] = spikenet2mex(p, act);
%       p -- the parameter structure filled in by spikenet1
%       act -- the activity array to overwrite with new activity.
%   returns:
%       uindex -- array of unit indexes (all in range 1..p.units)
%       tstamp -- array of matching timestamps.
%
% Runs the network for 10 seconds returns timestamps of firing units, and
% fills in the act array of unit activities.
 *
 * Version2 adds LFP calculations for each column.
 * Version7 changes conditioning to paired pulse instead of spike triggered.
 * Version8 adds inhibitory units to the paired pulse stimulation and allows
 *  for grading a units bias input.
 * Version9 adds weight decay and returns to the stimulus triggered paradigm.
 * Version10 allows a limit on how small a weight change can be made.
 * Version11 adds weight dependent scaling of wieght changes.
 * Version14 adds support for correlated bias input.
 * Version15 experiments with conditioning stimulations only on excitatory units
 * Version18 adds input rate modulation by a function to control column bursting
 * Version20 paired pulse conditioning.
 * Version25 back to spike triggered conditioning A(rec_unit)-> ColB
 * Version26 Synchrony handled as a percentage of bias spikes governed by bias_synchrony parameter
 * Version27 Timestep updates at t+1 rather than updating current activities and weights
 * Version30 Adjusts when STDP rules are applied, adds queue for spikes, adds correlated biases, removes synchrony percentage.
 * Version31 paired pulse conditioning.
 * Version32 Conditioning type is now a parameter.  Cycle triggered modulations added.
 * Version33 25% correlated biases / 75% uncorrelated.  Tetanic and Paired Pulse stims occur in between modulation episodes.
 * Version34 Fixed a bug in the axonslow[] update when training was turned off.  This bug allowed a much larger
 *  value for train_neg_factor due to accumulation of axonslow[unit] on the graph building trials which would give a big increase
 *  to weights at the beginning of the next training trial.  Maybe this worked a bit like simulated annealing?  Might have to
 *  try this out again, but it does explain some of the odd behaviors in the networks where weight seemed to be larger than expected.
 * Version35 saves an STDP time delta distribution for each weight A->B, A->C
 * Version36 Conditioning_Type 7 to trigger off of LFP when it rises through 0 after falling below p.lfp_detect_level.
 *  LFP is now the sum of all psps occurring in a column.   LFP[6] stores the filtered LFP for Col B
 *  Filters the 10kHz LFP with a 1st order Butterworth bandpass from [15Hz to 25Hz] 
 * Version37 Allows playback of previously recorded stimulation times.  This currently only makes sense for spike triggered
 *  and cycle triggered conditioning.
 * Version38 correlated bias inputs can now be modulated along with the uncorrelated bias inputs.
 * Version40 code cleanup. Conditioning type 3 replaced with the old conditioning type 7.
 * Version41 Activity Dependent Stimulation (conditioning type 7) changed to synchronus spikes in Col A causing delayed stim on Col B
 * Version42 Gamma Triggered Stimulation (conditioning type 8).
 * Version43 Spikes to output units have a second delay queue to increase delay time to upto 12 ms.
 * Version44 for Motor Output conditioning units in Column A create a local EMG in lfp[7].  Strength of the Ao1..AoN connections ramp.
 *   lfp[6] is a filtered 100-2500 Hz version to create action potentials that look more like motor unit action potentials.
 * Version 45 EMGA,B,C are stored in lfp[7], lfp[8], lfp[9].
 *   p.test_pulse_train is the number of pulses in the testing train.
 *   9/20/2019 gamma_filter_b gamma_filter_a to allow any band for the gamma triggered stimulation
 * Versoin 50 Allows timestep widths from 0.1 ms down to 0.01 ms in 0.1 ms steps.
 *   smaller timestep widths allow better approximations of the PSP and STDP functions.
 *   Conduction delay of spikes is now handled by a large circular buffer allowing a maximum 12ms delay at the 0.01 ms timestep
********************************************************************************/

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include <stdio.h>
#include <string.h>
#include <io.h>
#include <fcntl.h>
#include <stdint.h>

/* Input Arguments */

#define	PM prhs[0]
#define	ACT prhs[1]
#define	LFP prhs[2]

/* Output Arguments */

#define	UINDEX plhs[0]
#define TSTAMP plhs[1]

/* Constants */

// At the smallest timestep of 0.01 milliseconds: 
//   10 second simulation time is 1 million timsteps,
#define MAXCLOCK 1000000

// Increase MAXSPIKES if simulation complains about spike storage overflow
#define MAXSPIKES 1000000

// Increase MAXUNITS to accomodate larger networks
#define MAXUNITS 1000

// Maximum number of field potentials (or like things) we can track.
#define MAXLFP 100

/* prototypes */

typedef char int8;
typedef short int16;
typedef unsigned short uint16;
typedef int int32;
typedef unsigned int uint32;

/* Functions for returning information from a structure */

uint32 getM(const mxArray *pm, const char *fieldName)
{ // Number of rows in a field
    mxClassID classID;
    mxArray *fieldArray = mxGetField(pm, 0, fieldName);
    if (fieldArray != NULL) {
        classID = mxGetClassID(fieldArray);
        return (mxGetM(fieldArray)); // Max rows
    }
    return 0;    
}

uint32 getN(const mxArray *pm, const char *fieldName)
{ // Number of columns in a field
    mxClassID classID;
    mxArray *fieldArray = mxGetField(pm, 0, fieldName);
    if (fieldArray != NULL) {
        classID = mxGetClassID(fieldArray);
        return (mxGetN(fieldArray)); // Max rows
    }
    return 0;    
}

uint32 getLen(const mxArray *pm, const char *fieldName)
{ // Total number of elements stored in a field
    mxClassID classID;
    mxArray *fieldArray = mxGetField(pm, 0, fieldName);
    if (fieldArray != NULL) {
        classID = mxGetClassID(fieldArray);
        return (mxGetM(fieldArray) * mxGetN(fieldArray)); // Max elemets
    }
    return 0;    
}

void *getPr(const mxArray *pm, const char *fieldName)
{ // Get pointer to the underlying array.  Take care to match the
  // data type as there is no type checking or type conversion here.
    mxArray *fieldArray = mxGetField(pm, 0, fieldName);
    return mxGetPr(fieldArray);
}

double getValue(const mxArray *pm, const char *fieldName, int index)
{ // Gets Matlab equivalent of p.fieldName(index + 1)
    mxClassID classID;
    mxArray *fieldArray;
    int16 *i16;
    uint16 *u16;
    int32 *i32;
    uint32 *u32;
    uint32 len;
    double *doub;
    float *sing;
    
    fieldArray = mxGetField(pm, 0, fieldName);
    if (fieldArray != NULL) {
        classID = mxGetClassID(fieldArray);
        len = mxGetM(fieldArray) * mxGetN(fieldArray); // Max elemets
        if ((index < 0) || (index > len)) {
            index = len;  // out of bounds always returns last element
        }
        switch (classID) {
            case mxINT16_CLASS: 
                i16 = (int16 *)mxGetPr(fieldArray);
                return i16[index];
            case mxUINT16_CLASS: 
                u16 = (uint16 *)mxGetPr(fieldArray);
                return u16[index];
            case mxINT32_CLASS: 
                i32= (int32 *)mxGetPr(fieldArray);
                return i32[index];
            case mxUINT32_CLASS: 
                u32 = (uint32 *)mxGetPr(fieldArray);
                return u32[index];
            case mxDOUBLE_CLASS:
                doub = (double *)mxGetPr(fieldArray);
                return doub[index];
            case mxSINGLE_CLASS:
                sing = (float *)mxGetPr(fieldArray);
                return sing[index];
            default:
                mexPrintf("Field type not supported: %s\n", fieldName);
                return 0;
        }
    }
    mexPrintf("Field not found: %s\n", fieldName);
    return 0;
}

/* Main function */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])    
{ 
    // Static variables for storage and to maintain the state of the network across calls
    static int16 ulist[MAXSPIKES];
    static int32 tlist[MAXSPIKES];
    static double slowexp[MAXUNITS];
    static double fastexp[MAXUNITS];

    static double axonslow[MAXUNITS];  // Connection training decay function modifies connection strength when unit fires (strengthing rule)
    static double axonfast[MAXUNITS];  // (we can store this as a unit property since all connections have the same delay, this must change if we allow different connection delays)
    static double axondw[MAXUNITS];    // axonslow - axonfast for current time step.
    static double dendslow[MAXUNITS];  // Decay function at the dendrites (weakening rule)
    static double dendfast[MAXUNITS];  // 
    static double denddw[MAXUNITS];    // dendslow - dendfast for current time step
    static int32 stim_replay_index;    // Next clocktick to use from replay list.
    static int32 sim_clock=0; // Clock ticks counted from simulation start

    // Bias buffer size must be a power of 2.  It is the largest number of timesteps we can set a bias spike to occur in the future.
    // The size must be at least 8 * correlated_bias_std_ms * timesteps_per_second / 1000
    #define BIAS_BUFFER_MASK (0xfff)
    static int16 bias_buffer[MAXUNITS][BIAS_BUFFER_MASK+1]; 
    
    // Delay buffer size must be a power of 2.  It must be able to handle the maximum conduction time at the 0.01 ms timestep.
    // For example, a DELAY_BUFFER_MASK of 0x3ff gives a maximum conduction time of (2^10 - 1) * 0.01 = 10.23 milliseconds
    #define DELAY_BUFFER_MASK (0x3ff)
    static int8 spike_delay_buffer[MAXUNITS][DELAY_BUFFER_MASK+1]; // Circular buffer holding delayed spikes from each unit
    
    static double slowlfp[MAXLFP];
    static double fastlfp[MAXLFP];
    static uint32 rt;   // Random number generator state
    static uint32 rx;   // Random number generator state
    static uint32 ry;   // Random number generator state
    static uint32 rz;   // Random number generator state
    static uint32 rw;   // Random number generator state

    // Keep track of filtered LFPB state
    static double lfpX0;
    static double lfpX1;
    static double lfpX2;
    static double lfpY0;
    static double lfpY1;
    static double lfpY2;
    int lfpstate = 0; // 0 looking for LFP below threshold, 1 looking for LFP rising above 0;
    
    uint32 random_seed;   // Random seed for random number generator.
    int axonal_delay;     // Global axonal delay in timesteps.
    int dendritic_delay;  // Global dendritic delay in timesteps.
    int connection_delay; // Unit to Unit connection delay = axonal_delay + dendslow.
    int output_delay;     // Separate delay for connections to output units.
    int spike_queue_now_index; // Spike queue index for current clock tick
    int spike_queue_connection_delay_index; // .. at connection_delay
    int spike_queue_output_delay_index;     // .. at output_delay      
    
    int itrial;             // Trial index (1 = do initialization)
    int ntrials_before_test;// Number of trials before testing.
	double *bias_strength;  // = zeros(4, 1, 'int16'); % 
	uint32 *bias_chance;    // = zeros(4, p.time_steps, 'uint32'); % probability(step) = chance(step)/(2^32-1)
    double *train_info;     // place to keep track of training information
    int rec_unit;           // Index of the Record unit for ICMS, -1 to disable spike triggered stimulation.
    int rec_unit_flag = 0;  // Set to a positive number of timesteps when rec_unit fires.  Negative values used for counting out stim_refractory period.
    int rec_unit_count = 0; // Number of pulses left in the stimulus for spike triggered conditioning.
    int conditioning_flag;  // 0 for conditioning off, 1 for on.
    int conditioning_type;  //  % 0=No conditioning, 1=spike triggered, 2=paired pulse, 3 = cycle triggered LFP-B, 4 = Tetanic stim on Col A, 5 = Tetanic stim on Col B, 6 = Tetanic stim on Col C
    int train_flag;         // 0 = no training, 1= training allowed.
    int stim_delay;         // timesteps between rec_unit threshold crossing and stimulation on stim target units.
    int stim_pulse_train;   // Number of pulses in the spike triggered or paired pulse stimulus.
    int stim_refractory;    // minimum timesteps between stimulations for spike triggered and LFPB cycle triggered conditioning.
    double stim_phase_sine; // Sin() of the cycle triggered stimulation phase.
    int conditioning_end_time;// Time step to stop conditioning.
    int modulation_end_time;  // Time step to stop cycle triggered condioning (1 second longer than others for historical reasons)
    int spike_type = 0;     // 1 = spike goes to to non-output units, 2 = spike goes to output units, 3 = spike goest to to both
    
    double stim_uV;         // Stimulation size for conditioning.
    double pair_uV;         // Size of the paired pulse simulation (the second pulse, this is only for paired pulse stimulation)
    double test_uV;         // Stimulation for LFP testing.
    double max_psp_value;   // Limit on weight size
    double min_psp_value;   // Limit on inhibitory weight size.
    double weight_dependence; // weight dependence value. 1 for linear dependence, < 1 for increased bowing.
    
    double psp_slow_decay;    // Exponential decays governing psp shape
    double psp_fast_decay;
    
    double train_pos_factor;   // Hebbian learning rule for positive ms = (tPost - tPre)
    double train_pos_fast_decay;    // fast decay for time dependent learning rule.
    double train_pos_slow_decay;    // Slow decay for time dependent learning rule.

    double train_neg_factor;   // Hebbian learning rule for negative ms = (tPost - tPre)
    double train_neg_fast_decay;    // fast decay for time dependent learning rule.
    double train_neg_slow_decay;    // Slow decay for time dependent learning rule.

    double *unit_threshold;   // = zeros(p.units, 1, 'double');  % Each unit has its own threshold
    uint16 *unit_bias_offset; // = zeros(p.units, 1, 'uint16');  % Index of unit's bias
    uint16 *unit_pre_count;   // = zeros(p.units, 1, 'uint16');  % Number of presynaptic units
	uint16 *unit_post_count;  // = zeros(p.units, 1, 'uint16');  % Number of postsynaptic units
	uint32 *unit_pre_offset;  // = zeros(p.units, 1, 'uint32');  % Zero based offset into weight_pre_sort
    uint32 *unit_post_offset; // = zeros(p.units, 1, 'uint32');  % Zero based offset into wieght_post_sort
    uint16 *unit_lfp_offset;  // = zeros(p.units, 1, 'uint16');  % Zero based offset into LFP array
    uint16 *unit_column_id;   // = zeros(p.units, 1, 'uint16'); % 0=Col A unit, 1 = B, 2 = C, output units have 16384 added.
    uint16 *normal_pdf_table; // Contains 2^16 entries of normally distributed delay values for correlated input biases.
    double lfp_detect_level;  // Used for conditioning type 7 to detect oscillations in the lfp
    double lfp_amplitude_detect;  // The amplitude for detecting the next lfpb cycle triggered stimulus.
    double *gamma_filter_b;       // [b, a] filter coefficients for gamma triggered stimulation
    double *gamma_filter_a;
    double *beta_filter_b;         // LFP_B filter for cycle triggered stimulation
    double *beta_filter_a;
    double *emg_filter_b;         // Column A motor filed filter for emg triggerd stimulation
    double *emg_filter_a;
    double *unit_output_field_strength; // Motor field output strength values for constructing a synthethic motor output field.
    
    uint16 *unit_stim_source;  // flag for each unit if it is used for a leading paired pulse.
    uint16 *unit_stim_target;  // flag for each unit if it is used for a trailing paired pulse.
    uint16 *stim_source_times; // Stimulation flags for each clock tick (leading paired pulse)
    uint16 *stim_target_times; // Stimulation flags for each clock tick (trailing paired pulse)
    uint16 *stim_test_times;   // Times where we apply a test stimulation to 1 = Col A, 2 = Col B, 3 = Col C
    uint16 *stim_output_times; // Stimulation output flags for each clock tick
    uint16 *trigger_times;     // Stimulation output flags for each clock tick
    int stim_pulse_isi;        // Timesteps between stimulus pulses in a pulse train.
    int32 *stim_replay_clockticks;   // Clock ticks for forced conditioning stims.  This array must be terminated with a negative value.
    bool stim_replay_active = false; // True if we are replaying stimulations from the replay list.
    bool stim_replay_flag = false;   // True if current clock tick is time for a replayed stimulation.
            
    uint16 *weight_pre_unit;  // = zeros(p.weights, 1, 'uint16');      % Index of presynaptic unit (1..p.units)
	uint16 *weight_post_unit; // = zeros(p.weights, 1, 'uint16');      % Index of postsynpatic unit
	double *weight_strength;  // = zeros(p.weights, 1, 'double');      % Current connection strengths
	double *weight_strength2;  // = zeros(p.weights, 1, 'double');     % Connection strengths for next time step
	uint16 *weight_training_rule; // = zeros(p.weights, 1, 'uint16');  % Learning rule type (0 = none);
	uint32 *weight_pre_sort;  // = zeros(p.weights, 1, 'uint32');  % Offsets to weights sorted by pre_unit
	uint32 *weight_post_sort; // = zeros(p.weights, 1, 'uint32');  % Offsets to weights sorted by post_unit
    uint32 *weights_A_to_B;   // List of indexices of weights connecting A->B
    uint32 nWeights_A_to_B;   // Length of weights_A_to_B[];
    double *sumWeights_A_to_B; // Saves sum of mean weight A->B at each time step.
    uint16 *weight_test_lesion;  // Flag for each weight, 1 if testing should be done weight set to 0.
    
    double *act;    // Pointer to current element of the activity matrix    
    double *lfp;    // Pointer to current element of the lfp matrix
    int32 clock;    // keep track of time in clock ticks.
    int nCols;      // Number of columns
    int nExcit;    // Number of excitatory units in each column.
    int nInhib;     // Number of inhibitory units in each column.
    int nOut;       // Number of output units in each column.
    int nColUnits;  // nExcit + nInhib + nOut.
    int nSpikes;    // Number of spikes saved so far.
    int nUnits;     // Number of units in activity matrix.
    int nSteps;     // Number of time steps in activity matrix.
    int nWeights;   // Number of weights.
    int nBiases;    // Number of different types of biases.
    int nLfp;       // Number of LFP arrays.
    int biasid;     // Zero based bias id for current unit.
    int16 biaspsp;  // Strength of bias for current unit.
    int bias_offset;// offset to bias chances for current timestep.
    int unit;       // Unit index of current unit.
    int i;          // General index variable.
    int errFlag=0;  // Error flag
    int ilfp;       // offset of lfp for a unit
    uint16 iCol;    // Current unit's column id.
    int iu;
    int first_stim;
    int last_iCol;   // Used to check for first unit in a column.
	double *pUINDEX;		/* Access to the return variables. */
	double *pTSTAMP;
 
	/* Require 1 input parameter and 2 output parameters */

	if ((nlhs != 2) ||(nrhs != 3)) {
		mexErrMsgTxt("Usage: [uindex, tstamp] = spikenet1mex(p, act)"); 
		return;
	}
    
	/* Get input arguments */

    act = mxGetPr(ACT);  // Pointer to current unit's current activity.
    lfp = mxGetPr(LFP);  // Pointer to current LFP's current activity.
    nLfp = mxGetM(LFP);   // Number of LFP sums to track.
    nUnits = mxGetM(ACT); // Number of units
    nSteps = mxGetN(ACT); // Number of time steps to do
    nCols = getValue(PM, "n_cols", 0);
    nExcit = getValue(PM, "n_excit", 0);
    nInhib = getValue(PM, "n_inhib", 0);
    nOut = getValue(PM, "n_out", 0);
    nColUnits = nExcit + nInhib + nOut;
    nWeights = getValue(PM, "weights", 0); // Number of weights
    axonal_delay = getValue(PM, "axonal_delay", 0);
    dendritic_delay = getValue(PM, "dendritic_delay", 0);
    connection_delay = axonal_delay + dendritic_delay;
    output_delay = getValue(PM, "output_delay", 0);
    itrial = getValue(PM, "itrial", 0);
    ntrials_before_test = getValue(PM, "ntrials_on", 0) + getValue(PM, "ntrials_off", 0);
    rec_unit = getValue(PM, "rec_unit", 0) - 1; // zero based index of recording unit.
    train_flag = getValue(PM, "train_on", 0);   // 1 if training allowed.
    conditioning_flag = getValue(PM, "conditioning_flag", 0); // 1 if conditioning is on, 0 if off.
    conditioning_type = getValue(PM, "conditioning_type", 0);
    psp_fast_decay = getValue(PM, "psp_fast_decay", 0);
    psp_slow_decay = getValue(PM, "psp_slow_decay", 0);
    stim_delay = getValue(PM, "stim_delay", 0);
    stim_pulse_train = getValue(PM, "stim_pulse_train", 0);
    if (stim_pulse_train < 1) {
        stim_pulse_train = 1;
    }
    stim_pulse_isi = getValue(PM, "stim_pulse_isi", 0);
    gamma_filter_b = (double *)getPr(PM, "gamma_filter_b");       // LFP_A filter for gamma triggered stimulation
    gamma_filter_a = (double *)getPr(PM, "gamma_filter_a");
    beta_filter_b = (double *)getPr(PM, "beta_filter_b");         // LFP_B filter for cycle triggered stimulation
    beta_filter_a = (double *)getPr(PM, "beta_filter_a");
    emg_filter_b = (double *)getPr(PM, "emg_filter_b");          // Column A motor filed filter for emg triggerd stimulation
    emg_filter_a = (double *)getPr(PM, "emg_filter_a");
    
    stim_refractory = -getValue(PM, "stim_refractory", 0);
    stim_phase_sine = getValue(PM, "stim_phase_sine", 0);
    conditioning_end_time = nSteps * getValue(PM, "conditioning_secs", 0) / 10;
    modulation_end_time = conditioning_end_time + (nSteps / 10);
    stim_uV = getValue(PM, "stim_uV", 0);
    pair_uV = getValue(PM, "pair_uV", 0);
    test_uV = getValue(PM, "test_uV", 0);
    max_psp_value = getValue(PM, "max_psp_value", 0);
    min_psp_value = -max_psp_value;
    weight_dependence = getValue(PM, "weight_dependence", 0);
    random_seed = (uint32)getValue(PM, "random_seed", 0);

	bias_strength = (double *)getPr(PM, "bias_strength"); 
    bias_chance = (uint32 *)getPr(PM, "bias_chance");
    nBiases = getM(PM, "bias_chance");
           
    train_info = (double *)getPr(PM, "train_info");    
    train_pos_factor = getValue(PM, "train_pos_factor", 0);
    train_pos_slow_decay = getValue(PM, "train_pos_slow_decay", 0);
    train_pos_fast_decay = getValue(PM, "train_pos_fast_decay", 0);
    train_neg_factor = getValue(PM, "train_neg_factor", 0);
    train_neg_slow_decay = getValue(PM, "train_neg_slow_decay", 0);
    train_neg_fast_decay = getValue(PM, "train_neg_fast_decay", 0);
  
	unit_threshold = (double *)getPr(PM, "unit_threshold");
	unit_bias_offset = (uint16 *)getPr(PM, "unit_bias_offset");
	unit_pre_count = (uint16 *)getPr(PM, "unit_pre_count");
	unit_post_count = (uint16 *)getPr(PM, "unit_post_count");
	unit_pre_offset = (uint32 *)getPr(PM, "unit_pre_offset");
    unit_post_offset = (uint32 *)getPr(PM, "unit_post_offset");
    unit_lfp_offset = (uint16 *)getPr(PM, "unit_lfp_offset");
    unit_column_id = (uint16 *)getPr(PM, "unit_column_id");
    normal_pdf_table = (uint16 *)getPr(PM, "normal_pdf_table");
    lfp_detect_level = getValue(PM, "lfp_detect_level", 0);
    unit_output_field_strength = (double *)getPr(PM, "unit_output_field_strength");
    
	weight_pre_unit = (uint16 *)getPr(PM, "weight_pre_unit");
	weight_post_unit = (uint16 *)getPr(PM, "weight_post_unit");
	weight_strength = (double *)getPr(PM, "weight_strength");
	weight_training_rule = (uint16 *)getPr(PM, "weight_training_rule");
	weight_pre_sort = (uint32 *)getPr(PM, "weight_pre_sort");
	weight_post_sort = (uint32 *)getPr(PM, "weight_post_sort");
    sumWeights_A_to_B = (double *)getPr(PM, "sumWeights_A_to_B");
    weights_A_to_B = (uint32 *)getPr(PM, "weights_A_to_B");
    nWeights_A_to_B = getM(PM, "weights_A_to_B");
	weight_test_lesion = (uint16 *)getPr(PM, "weight_test_lesion");
            
    // Only use one of the following two lines.
	//weight_strength2 = (double *)getPr(PM, "weight_strength2");  // Update weights at time t+1.  Slower but equation accurate.
	weight_strength2 = weight_strength;  // Optimization to update weights in place at time t.  Faster and nearly the same results.
    
    // Flags for paired pulse stimulation
    unit_stim_source = (uint16 *)getPr(PM, "unit_stim_source");
    unit_stim_target = (uint16 *)getPr(PM, "unit_stim_target");
    stim_source_times = (uint16 *)getPr(PM, "stim_source_times");
    stim_target_times = (uint16 *)getPr(PM, "stim_target_times");
    stim_test_times = (uint16 *)getPr(PM, "stim_test_times");
    stim_output_times = (uint16 *)getPr(PM, "stim_output_times");
    trigger_times = (uint16 *)getPr(PM, "trigger_times");
    stim_replay_clockticks = (int32 *)getPr(PM, "stim_replay_clockticks");
    
    //mexPrintf("Units: %g, Steps: %g, Weights: %g, delay: %g, biases: %g\n", (double)nUnits, (double)nSteps, (double)nWeights, (double)axonal_delay, (double)nBiases);
    //mexPrintf("stimdelay: %g\n", (double)stim_delay);
    
    /* Initialize activity on first sweep */
    
    if (itrial == 1) {
        int i;
        sim_clock = 0;
        for (unit=0; unit < MAXUNITS; unit++) {
            slowexp[unit] = 0;
            fastexp[unit] = 0;
            dendslow[unit] = 0;
            dendfast[unit] = 0;
            axonslow[unit] = 0;
            axonfast[unit] = 0;
            for (i=0; i <= BIAS_BUFFER_MASK; i++) {
                bias_buffer[unit][i] = 0;
            }
            for (i=0; i <= DELAY_BUFFER_MASK; i++) {
                spike_delay_buffer[unit][i] = 0;
            }
        }
        
        rw = 1; rx = 10; ry = random_seed; rz = 1000; // Init RNG state.
        
        lfpX0 = 0;  lfpX1 = 0;  lfpX2 = 0;
        lfpY0 = 0;  lfpY1 = 0;  lfpY2 = 0;
        lfpstate = 0; // 0 looking for LFP below threshold
        stim_replay_index = 0;    // Next clocktick to use from replay list.
    }
        
   /* Run network */
 
    nSpikes = 0; // Clear spike list
    bias_offset = 0;  // Starting at beginning of bias_chance matrix.
    stim_replay_active = (stim_replay_clockticks[stim_replay_index] > 0);
    for (clock = 0; clock < nSteps; clock++, sim_clock++, bias_offset += nBiases)
    {
        if (stim_replay_active) { // Check if we are replaying stimulations
            stim_replay_flag = (stim_replay_clockticks[stim_replay_index] == (itrial - 1) * nSteps + clock);
        }
        
        // Precalculate indexes into spike_delay_buffer
        spike_queue_now_index = sim_clock & DELAY_BUFFER_MASK; // Place to put spikes occuring at this clock tick
        spike_queue_connection_delay_index = (sim_clock - connection_delay) & DELAY_BUFFER_MASK;  // Place where spikes occurred at normal connection delay
        spike_queue_output_delay_index = (sim_clock - output_delay) & DELAY_BUFFER_MASK;  // Place where spikes occurred at the output connection delay      
        
        // Update unit activity for this time step
        for (unit = 0; unit < nUnits; unit++)
        {
            act[unit] = slowexp[unit] - fastexp[unit];   // Activity of unit at this time step.
            slowexp[unit] *= psp_slow_decay;  // Slow decay potential for next time step
            fastexp[unit] *= psp_fast_decay;  // Fast decay potential for next time step
            if (train_flag > 0) {
                axondw[unit] = axonslow[unit] - axonfast[unit]; // Current weight change for strenthening rule.
                axonslow[unit] *= train_pos_slow_decay;  // Slow decay of connection firing training potential
                axonfast[unit] *= train_pos_fast_decay;  // Fast decay of connection firing training potential
                denddw[unit] = dendslow[unit] - dendfast[unit]; // Current weight change for weaking rule.
                dendslow[unit] *= train_neg_slow_decay;  // Slow decay of unit firing training potential
                dendfast[unit] *= train_neg_fast_decay;  // Fast decay of unit firing training potential
            }
        }
        
        // Update output for all LFP arrays for this timestep.
        for (i = 0; i < nLfp; i++) {
            *lfp++ = slowlfp[i] - fastlfp[i]; // resulting LFP at this timestep
            slowlfp[i] *= psp_slow_decay;  // Slow decay for next timestep
            fastlfp[i] *= psp_fast_decay;  // Fast decay for next timestep
        }
        
        // For motor output conditioning, filter the ColA motor field with a 100-2500 Hz filter
        if (conditioning_type == 7) {
            lfpX2 = lfpX1;
            lfpX1 = lfpX0;
            lfpY2 = lfpY1;
            lfpY1 = lfpY0;
            lfpX0 = slowlfp[7] - fastlfp[7]; // Get current EMGA
            
            lfpY0 = (lfpX0 - lfpX2) * emg_filter_b[0] - lfpY1 * emg_filter_a[1] - lfpY2 * emg_filter_a[2];
            slowlfp[6] = lfpY0;  // Save filtered Column A output unit LFP
            fastlfp[6] = 0;
            
            if (lfp_detect_level < 0) {
                // Handle falling motor unit field
                if (lfpY0 < lfp_detect_level) {  // Wait for filtered LFPA to go below detection level
                    if (lfpstate == 0) {
                        lfpstate = 1;
                        if ((rec_unit_flag <= stim_refractory) && !stim_replay_active) { // Check for refractory period
                            rec_unit_flag = stim_delay + 1;  // Stimulate target units after delay
                            rec_unit_count = stim_pulse_train; // Number of stimulus pulses to deliver
                            trigger_times[clock] = 1;         // Remember time of trigger detection
                        }
                    }
                } else if (lfpY0 > 0) {
                    lfpstate = 0; // Wait for filtered lfp to go above 0 before allowing another trigger
                }
            } else {
                // Handle rising motor unit field
                if (lfpY0 > lfp_detect_level) {  // Wait for filtered LFPA to go above detection level
                    if (lfpstate == 0) {
                        lfpstate = 1;
                        if ((rec_unit_flag <= stim_refractory) && !stim_replay_active) { // Check for refractory period
                            rec_unit_flag = stim_delay + 1;  // Stimulate target units after delay
                            rec_unit_count = stim_pulse_train; // Number of stimulus pulses to deliver
                            trigger_times[clock] = 1;         // Remember time of trigger detection
                        }
                    }
                } else if (lfpY0 < 0) {
                    lfpstate = 0; // Wait for filtered lfp to go below 0 before allowing another trigger
                }
            }
        } // if (conditioning_type == 7)
        
        // Update calculaton for the Gamma filtered LFPA when conditioning type 8 is used,
        if (conditioning_type == 8) {
            lfpX2 = lfpX1;
            lfpX1 = lfpX0;
            lfpY2 = lfpY1;
            lfpY1 = lfpY0;
            lfpX0 = (slowlfp[0] - fastlfp[0]) + (slowlfp[3] - fastlfp[3]);
            
            // Mid gamma 50-80 Hz: [b, a] = butter(1, [50 80] / fs/2); num2str([b(1) a(2) a(3)], 20)
            // Low gamma 30-50 Hz: [b, a] = butter(1, [30 50] / fs/2); num2str([b(1) a(2) a(3)], 20)
            // High gamma 80-100 Hz: [b, a] = butter(1, [80 100] / fs/2); num2str([b(1) a(2) a(3)], 20)
            lfpY0 = (lfpX0 - lfpX2) * gamma_filter_b[0] - lfpY1 * gamma_filter_a[1] - lfpY2 * gamma_filter_a[2];
            slowlfp[6] = lfpY0;  // Save filtered LFPA
            fastlfp[6] = 0;
            
            if (lfp_detect_level < 0) {
                // Handle falling LFPA (for all negative values of lfp_detect_level)
                if (lfpY0 < lfp_detect_level) {  // Wait for filtered LFPA to go below detection level
                    if (lfpstate == 0) {
                        lfpstate = 1;
                        if ((rec_unit_flag <= stim_refractory) && !stim_replay_active) { // Check for refractory period
                            rec_unit_flag = stim_delay + 1;  // Stimulate target units after delay
                            rec_unit_count = stim_pulse_train; // Number of stimulus pulses to deliver
                            trigger_times[clock] = 1;         // Remember time of trigger detection
                        }
                    }
                } else if (lfpY0 > 0) {
                    lfpstate = 0; // Wait for filtered lfp to go above 0 before allowing another trigger
                }
            } else {
                // Handle rising LFPA (for all non-negative values of lfp_detect_level)
                if (lfpY0 > lfp_detect_level) {  // Wait for filtered LFPA to go above detection level
                    if (lfpstate == 0) {
                        lfpstate = 1;
                        if ((rec_unit_flag <= stim_refractory) && !stim_replay_active) { // Check for refractory period
                            rec_unit_flag = stim_delay + 1;  // Stimulate target units after delay
                            rec_unit_count = stim_pulse_train; // Number of stimulus pulses to deliver
                            trigger_times[clock] = 1;         // Remember time of trigger detection
                        }
                    }
                } else if (lfpY0 < 0) {
                    lfpstate = 0; // Wait for filtered lfp to go below 0 before allowing another trigger
                }
            }
        }
               
        // Update calculation for the cycle triggered band filtered LFPB when conditioning type 3 is used.
        if (conditioning_type == 3) {
            lfpX2 = lfpX1;
            lfpX1 = lfpX0;
            lfpY2 = lfpY1;
            lfpY1 = lfpY0;
            lfpX0 = (slowlfp[1] - fastlfp[1]) + (slowlfp[4] - fastlfp[4]);
            lfpY0 = (lfpX0 - lfpX2) * beta_filter_b[0] - lfpY1 * beta_filter_a[1] - lfpY2 * beta_filter_a[2];
            slowlfp[6] = lfpY0;  // Save filtered LFPB
            fastlfp[6] = 0;
            
            if (clock < modulation_end_time) {
                if (lfp_detect_level < 0) {
                    // Handle stim phases from -90 to less than 90
                    if (lfpstate != 1) {
                        if (lfpY0 < lfp_detect_level) {  // Wait for filtered LFPB to go below detection level
                            if (lfpY1 < lfpY0) {  // Wait for peak (hopefully the filter removes most of the jitter here)
                                lfpstate = 1;
                                lfp_amplitude_detect = stim_phase_sine * -lfpY0;
                            }
                        }
                    }
                    if (lfpstate != 0) {
                        if (lfpY0 >= lfp_amplitude_detect) {  // Wait for LFPB to cross our detection threshold
                            lfpstate = 0;
                            if ((rec_unit_flag <= stim_refractory) && !stim_replay_active) { // Check for refractory period
                                rec_unit_flag = 1;  // Stimulate target units on this timestep
                                rec_unit_count = stim_pulse_train; // Number of stimulus pulses to deliver
                                trigger_times[clock] = 1;         // Remember time of trigger detection
                            }
                        } else if (lfpY0 > 0) {
                            lfpstate = 2;   // Start looking at lfp_detect_level again
                        }
                    }
                } else {
                    // Handle stim phases from 90 to less than 270
                    if (lfpstate != 1) {
                        if (lfpY0 > lfp_detect_level) {  // Wait for filtered LFPB to go above detection level
                            if (lfpY1 > lfpY0) {  // Wait for peak
                                lfpstate = 1;
                                lfp_amplitude_detect = stim_phase_sine * lfpY0;
                            }
                        }
                    }
                    if (lfpstate != 0) {
                        if (lfpY0 <= lfp_amplitude_detect) {          // Then wait for LFPB to cross 0
                            lfpstate = 0;
                            if ((rec_unit_flag <= stim_refractory) && !stim_replay_active) { // Check for refractory period
                                rec_unit_flag = 1;  // Stimulate target units on this timestep
                                rec_unit_count = stim_pulse_train; // Number of stimulus pulses to deliver
                                trigger_times[clock] = 1;         // Remember time of trigger detection
                            }
                        } else if (lfpY0 < 0) {
                            lfpstate = 2;  // Start looking at lfp_detect_level again
                        }
                    }
                }
            }
        }
        
        // Update weights from last time step.
        // This slow things down a bit so there is an optimization that updates
        // weights in-place at time t instead of at time t+1.  While this
        // isn't an exact match to the equations, it yeilds nearly the
        // same results in a little less time.
        if (weight_strength != weight_strength2) {
            double *p1 = weight_strength;
            double *p2 = weight_strength2;
            int icount = nWeights;
            do {
                *p1++ = *p2++;
            } while (--icount);
        }
        
        // Handle bias inputs to each unit and check for unit spiking
        first_stim = 0;
        last_iCol = 100;  // Last column index used to check when unit group switches
        for (unit = 0; unit < nUnits; unit++) {
            biasid = unit_bias_offset[unit];
            ilfp = unit_lfp_offset[unit];    // Target unit's LFP index.
            iCol = unit_column_id[unit];
            
            // Check for correlated input biases.
            
            if (last_iCol != iCol) {  // Column index switched ...
                last_iCol = iCol;     // ... create correlated bias spikes for this column, but 
                if (iCol < 100) {     // ... don't apply this to output units (which are marked with large column ids)
                    rt = rx ^ (rx << 11);       // Update RNG state
                    rx = ry; ry = rz; rz = rw;
                    rw = rw ^ (rw >> 19) ^ rt ^ (rt >> 8);
                    if (rw < bias_chance[biasid + bias_offset]) { // For all units in this column                    
                        int destunit = unit;
                        int bin;
                        for (destunit = unit; unit_column_id[destunit] == iCol; destunit++) {
                            rt = rx ^ (rx << 11);           // Update RNG state
                            rx = ry; ry = rz; rz = rw;
                            rw = rw ^ (rw >> 19) ^ rt ^ (rt >> 8);
                            bin = (sim_clock + normal_pdf_table[rw & 0xffff]) & BIAS_BUFFER_MASK;
                            bias_buffer[destunit][bin] += 1;
                        }
                        //mexPrintf("unit %d, dest %d, col %d %d, biasid %d\n", unit, destunit, iCol, unit_column_id[destunit], biasid]);                    
                    }
                }
            }
            
            biaspsp = bias_buffer[unit][sim_clock & BIAS_BUFFER_MASK];
            if (biaspsp != 0) {
                // Correlated bias spike(s) occurs on this unit at this timestep.
                train_info[3] += biaspsp; // Count a correlated bias spike.
                bias_buffer[unit][sim_clock & BIAS_BUFFER_MASK] = 0; // Clear circular buffer entry.
            }
            
            // Check if it is time for an uncorrelated bias input spike on this unit.
            // Uncorrelated chance is must be stored 1 index after the correlated chance.
            
            rt = rx ^ (rx << 11);               // Update RNG state
            rx = ry; ry = rz; rz = rw;
            rw = rw ^ (rw >> 19) ^ rt ^ (rt >> 8);
            if (rw < bias_chance[biasid + 1 + bias_offset]) {
                biaspsp++;
            }
            
            // Add in all bias psps to the unit potential
            
            if (biaspsp > 0) {
                double biasval = biaspsp * bias_strength[biasid];
                slowexp[unit] += biasval;  // Update slow decay potential with bias
                fastexp[unit] += biasval;  // Update fast decay potential with bias
                slowlfp[ilfp] += biasval;  // Handle the LFP calculation separately.
                fastlfp[ilfp] += biasval;
            }
           
            // Check if it is time to do a paired pulse, tetanic, activity triggered, or spike-triggered stimulation.
            
            if (conditioning_flag != 0) {
                if (unit_stim_source[unit] & stim_source_times[clock]) {
                    slowexp[unit] += stim_uV;  // Do initial paired pulse on the source units.                      
                    rec_unit_flag = 0;         // Reset refractory period
                    if (first_stim == 0) {
                        first_stim = 1;        // Count number of stimuli only on first unit encountered for each clock tick.
                        stim_output_times[clock] = 1;
                        train_info[2]++;
                        trigger_times[clock] = 1;   // Remember time of trigger detection
                    }
                } else if (unit_stim_target[unit] & stim_target_times[clock]) {
                    slowexp[unit] += pair_uV;  // Do delayed paired pulse on target units, or test pulses for lfp cyclic conditioning.
                    rec_unit_flag = 0;         // Reset refractory period
                } else if (((rec_unit_flag == 1) || stim_replay_flag) && (unit_stim_target[unit] != 0)) {
                    // Spike triggered, target paired pulse, or replayed stimulations are done here.
                    if (clock < conditioning_end_time) {                  
                        slowexp[unit] += stim_uV;  // Update slow decay potential with stimulus
                        if (first_stim == 0) {
                            first_stim = 1;        // Count number of stimuli only on first unit encountered for each clock tick.
                            stim_output_times[clock] = 1;
                            train_info[2]++;
                        }
                    }
                }
            } else {
                // Even when conditioning is off, keep track of trigger times and refractory period.
                if (unit_stim_source[unit] & stim_source_times[clock]) {
                    rec_unit_flag = 0;         // Reset refractory period
                    if (first_stim == 0) {
                        first_stim = 1;        // Count number of stimuli only on first unit encountered for each clock tick.
                        trigger_times[clock] = 1;   // Remember time of trigger detection
                    }
                } else if (unit_stim_target[unit] & stim_target_times[clock]) {
                    rec_unit_flag = 0;         // Reset refractory period
                } else if (((rec_unit_flag == 1) || stim_replay_flag) && (unit_stim_target[unit] != 0)) {
                    // Spike triggered, activity dependent, target paired pulse, or replayed stimulations are done here.
                    if (first_stim == 0) {
                        first_stim = 1;        // Count number of stimuli only on first unit encountered for each clock tick.
                    }
                }
            }
            
            if (train_flag == 0) {
                // Stimulus test pulses at specfic times during testing sections.
                if (stim_test_times[clock] == unit_column_id[unit] + 1) { // Check for test stim time on this unit's column.
                    slowexp[unit] += test_uV;  // Update slow decay potential with test stim on this unit.
                    //mexPrintf("Stim %d on col %d\n", clock, unit_column_id[unit] + 1);
                }
            }
            
            // Check for unit activity threshold crossings and connection firings
                      
            spike_delay_buffer[unit][spike_queue_now_index] = 0; // Clear place to hold spike occurrence
            
            if (act[unit] >= unit_threshold[unit]) { // Threshold crossed
                // Place spike in queue.  Clear unit potential.  Record spike time.
                spike_delay_buffer[unit][spike_queue_now_index] = 1;   // Start spike
                slowexp[unit] = 0;  // Reset slow and fast exponentials
                fastexp[unit] = 0;

                // For spike-triggered stimulation.
                // Fake stimulation on the target units in a certain number of timesteps
                if ((unit == rec_unit) && (rec_unit_flag <= stim_refractory) && !stim_replay_active) {
                    if (conditioning_flag != 0) {
                        rec_unit_flag = stim_delay + 1;    // Deliver a conditioning stim at a delay
                        rec_unit_count = stim_pulse_train; // Number of stimulus pulses to deliver
                    }
                    trigger_times[clock] = 1;         // Remember time of trigger detection
                }
                
                if (nSpikes < MAXSPIKES) {    // Store unit firing time
                    tlist[nSpikes] = clock;
                    ulist[nSpikes++] = unit;
                } else {
                    if (errFlag == 0) {
                        mexPrintf("Spike storage overflow at timestep %g\n", (double)clock);
                    }
                    errFlag = 1;
                }
                                
                if (train_flag > 0) {
                    // Handle training positive rule for connections to this unit that fired previously.
                    // Run through list of weights having the current unit as their post-synaptic unit.
                    int iwend;
                    int iw = unit_post_offset[unit];  // Offset into post_sort list
                    for (iwend = iw + unit_post_count[unit]; iw < iwend; iw++) {
                        int windex = weight_post_sort[iw];
                        
                        // sanity check
                        //if (weight_post_unit[windex]-1 != unit) {
                        //    mexPrintf("expecting unit %d != %d\n", unit, weight_pre_unit[windex]-1);
                        //}
                        
                        if (weight_training_rule[windex] > 0) {
                            double dw;
                            double weight = weight_strength[windex]; // Current weight strength.
                            iu = weight_pre_unit[windex] - 1;    // Convert to zero based array index.
                            dw = axondw[iu]; // Weight change for strengthening rule.
                            if (weight > 0) {
                                // Excitatory connection, apply weight change -> max_weight
                                dw *= pow(1 - (weight / max_psp_value), weight_dependence);
                                weight_strength2[windex] += dw;  // Strengthen excite conn when post spike follows.
                                if (weight_strength2[windex] > max_psp_value) {
                                    weight_strength2[windex] = max_psp_value;
                                }
                                //train_info[0] += dw;
                                train_info[0] += 1;
                            } else if (weight < 0) {
                                // Inhibitory connection, apply weight change -> min_weight.
                                dw *= pow(1 - (-weight / max_psp_value), weight_dependence);
                                weight_strength2[windex] -= dw;  // Strengthen inhib conn when post spike follows.
                                if (weight_strength2[windex] < min_psp_value) {
                                    weight_strength2[windex] = min_psp_value;
                                }
                                //train_info[1] += dw;
                                train_info[0] += 1;
                            }
                        }
                    }
                    
                    // sanity check
                    //if (weight_post_unit[weight_post_sort[iw + unit_post_count[unit]]]-1 == unit) {
                    //    mexPrintf("check != %d, %d\n", unit, weight_post_unit[weight_post_sort[iw + unit_post_count[unit]]]-1);
                    //}
                    
                    // Update training function for the negative rule (for when connections to this unit fire after this time)
                    dendslow[unit] += train_neg_factor;
                    dendfast[unit] += train_neg_factor;
                                        
                } // end if (train_flag

                // Outputs calculate a synthetic EMG
                if ((iCol >= 16384)) { // Output unit column is tagged with +16384
                    double mfs = unit_output_field_strength[unit]; // EMG contribution
                    int ind = iCol - 16384 + 7; // EMG A at index 7, EMG B at 8, EMG C at 9
                    slowlfp[ind] = slowlfp[ind] + mfs;  // This will normally be band pass filtered before use,
                    fastlfp[ind] = fastlfp[ind] + mfs;  // so it will end up lower amplitude in averages and such.         
                }
                
            } // end if (act[unit]
            
            // Check if this unit spiked at conection_delay or output_connection delay timesteps ago.
            spike_type = 0;
            if (spike_delay_buffer[unit][spike_queue_connection_delay_index]) {
                spike_type = 1;
            }
            if (spike_delay_buffer[unit][spike_queue_output_delay_index]) {
                spike_type |= 2;
            }
                            
            if (spike_type) {
                // Deliver psps from current unit to post synaptic units at this time.
                int iwend;
                int iw = unit_pre_offset[unit];  // Offset into pre_sort list

                for (iwend = iw + unit_pre_count[unit]; iw < iwend; iw++) {
                    int windex = weight_pre_sort[iw];
                    double weight = weight_strength[windex];
                    iu = weight_post_unit[windex] - 1; // Convert to zero based array index.
                    ilfp = unit_lfp_offset[iu];        // Post synaptic unit's LFP index. 0 ColA, 1 ColB, 2 ColC, 9 ColA Output, 10 ColB Output, 11 ColC Output
                    
                    // sanity check
                    // if (weight_pre_unit[windex]-1 != unit) {
                    //     mexPrintf("expecting pre unit %d != %d\n", unit, weight_pre_unit[windex]-1);
                    // }
                    
                    if ( ((ilfp <= 3) && ((spike_type & 1) != 0)) || ((ilfp >= 9) && ((spike_type & 2) != 0)) ) {

                       // Disallow PSPs for marked weights when plasticity is off (i.e. during testing).
                        if ((train_flag > 0) || (weight_test_lesion[windex] == 0)) {
                            // Psp to connecting unit 
                            slowexp[iu] += weight;
                            fastexp[iu] += weight;
                            //mexPrintf("Post unit %d, train rule %d\n", iu, weight_training_rule[windex]);

                            // Handle LFP calculation
                            if (weight > 0) {
                                slowlfp[ilfp] += weight;
                                fastlfp[ilfp] += weight;
                            } else { // Separate inhibitory contribution to LFP
                                slowlfp[ilfp + 3] += weight;
                                fastlfp[ilfp + 3] += weight;
                            }
                        }
                        
                        if (train_flag > 0) {
                            if (weight_training_rule[windex] > 0) {
                                double dw = denddw[iu]; // Weight change for weaking rule
                                //mexPrintf("unit %g, weight %g\n", (double)unit, (double)weight, (double)weight);
                                if (weight > 0) {
                                    // Excitatory connection, apply weight change -> zero.
                                    dw *= pow((weight / max_psp_value), weight_dependence);
                                    weight_strength2[windex] -= dw;
                                    if (weight_strength2[windex] < 1) {
                                        weight_strength2[windex] = 1; // Excitatory weight must stay >= 1
                                    }
                                    //train_info[1] += dw;
                                    train_info[1] += 1;
                                } else if (weight < 0) {
                                    // Inhibitory connection, apply weight change -> zero
                                    dw *= pow((-weight / max_psp_value), weight_dependence);
                                    weight_strength2[windex] += dw;  // Weaken excite conn when pre spike follows.
                                    if (weight_strength2[windex] > -1) {
                                        weight_strength2[windex] = -1; // Inhibitory weight must stay <= -1
                                    }
                                    //train_info[0] += dw;
                                    train_info[1] += 1;
                                }
                            } // end if (weight_training_rule
                        } // if (train_flag
                    } // if (ilfp and spike_type are a match
                } // end for (iwend
                
                // sanity check
               // if (weight_pre_unit[weight_pre_sort[iw + unit_pre_count[unit]]]-1 == unit) {
               //     mexPrintf("check != %d, %d\n", unit, weight_pre_unit[weight_pre_sort[iw + unit_pre_count[unit]]]-1);
               // }
                
                // Update training function for the positive rule (for when connections to this unit fire before this time)
                if ((train_flag > 0) && ((spike_type & 1) != 0)) {
                    axonslow[unit] += train_pos_factor;
                    axonfast[unit] += train_pos_factor;
                }
            } // end if (queue
            
        } // end for (unit 
        
        // Point to next activity time step        
        *act += nUnits;
        
        // Countdown timer for stimulation potentials from the recording unit.
        rec_unit_flag--;  // Decrement rec_unit_flag, negative values used for counting out stim_refractory period
        // Handle aditional pulses in spike triggered stimulus train.
        if ((rec_unit_flag == 0) && (rec_unit_count > 0)) {
            rec_unit_count--;
            if (rec_unit_count > 0) {
                rec_unit_flag = stim_pulse_isi; // Next pulse occurs at interstimulus interval
            }
        }
        
        if (stim_replay_flag) {
            stim_replay_index++; // Index for next stim replay clocktick.
        }
        
        // Find the mean of all weights connecting Column A to Column B.
        if (conditioning_type > 0) { // Possibly don't do this for all conditioning types?
            int iw;
            double sum = 0;
            for (iw=0; iw<nWeights_A_to_B; iw++) {
                sum += weight_strength[weights_A_to_B[iw]];
            }
            sumWeights_A_to_B[clock] += sum / nWeights_A_to_B;
        }
        
    } // end for (clock
    
   /* Allocate memory for return values. */

	UINDEX = mxCreateDoubleMatrix(nSpikes, 1, mxREAL);
	TSTAMP = mxCreateDoubleMatrix(nSpikes, 1, mxREAL);

	if (!UINDEX || !TSTAMP) {
		if (UINDEX) mxFree(UINDEX); UINDEX = NULL;
		if (TSTAMP) mxFree(TSTAMP); TSTAMP = NULL;
		mexErrMsgTxt("Error: out of memory in spikenetmex().");
		return;
	}

	pUINDEX = mxGetPr(UINDEX);
	pTSTAMP = mxGetPr(TSTAMP);
    for (i=0; i<nSpikes; i++) {
        *pUINDEX++ = ulist[i] + 1;
        *pTSTAMP++ = tlist[i] + 1;
    }

	return;
}
