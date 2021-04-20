% SPIKING-NETWORK-SIMULATIONS Based on Shupe & Fetz (2021) eNeuro paper.
%  See Larry Shupe's repo for original implementation:
%  >   https://github.com/lshupe/spikenet
%  This repo is forked from commit:
%  >  cfb29cefa80f61bc8bde97ec279300e2a3e50a69
%  The modifications in this repository are for three purposes:
%     1. Provide formatted documentation as I read/interpret their work.
%     2. Provide an object-oriented API to make modifications/use easier.
%     3. Reproduce their experimental results and see how robust they are
%        to specific parametric changes, such as the probability of
%        connection for synapses capable of generating specific kinds of
%        post-synaptic currents.
%  About Mex files:
%     Mex files are basically a way to use pointers in Matlab. 
%     This is useful to run certain things (such as simulations) more 
%     efficiently. For example, a lot of the methods in the mex file of 
%     this project involve using bit-wise operations to approximate math 
%     operations to allow the simulator to run quickly. 
%     Errors involving spikenet50mex.mexw64 mean you need to check
%     spikenet50mex.c. The mex file has special code specifically for
%     handling expected fields of a specific struct. So a probable cause of
%     error is differences in field naming conventions between field names
%     of variables in the Matlab code vs. what the c code is expecting in
%     terms of the struct field names. In `SpikeNet`, there is a
%     `get_mex_struct` method that should export the parameters into the
%     correct struct naming convention.
%
% Classes
%   SimulationData        - A class for saving output data from simulations.
%   SpikeNet              - Object-oriented implementation of spikenet50 (Larry Shupe)
%
% Functions
%  spikenet50_spike_trig  - Runs an integerate and fire neural network model.
%
% Mex
%  spikenet50mex.mexw64   - Compiled for Win10 64-bit, Matlab R2020b.
%  spikenet50mex.c        - Source file for compiled mex file.
%
% Packages
%  +f                     - Package containing functions with graphics utilities
%
% Folders
%  original               - Larry's original functions should be kept here.
%
% Scripts
%   main                  - This script is the main way that the classes are initialized / used
