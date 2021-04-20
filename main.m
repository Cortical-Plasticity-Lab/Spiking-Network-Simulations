%MAIN This script is the main way that the classes are initialized / used
clear; clc;

% CHANGE THESE PARAMETERS TO YOUR NEEDS %
LOCAL_DATA_SAVE_LOCATION = 'D:/data/spikenet'; % All outputs saved here!

% Any spikenet original functions of interest should be at the top-level of
% the repository, and would be called from here.
spikenet50_spike_trig(LOCAL_DATA_SAVE_LOCATION); 

% % % % WIP: This part does not work correctly (yet) % % % % %
% spikenet = SpikeNet(); % CONTROL "default" spikenet50
% 
% % Set up ADS spikenet (conditioning_type == 1)
% spikenet_ads = copy(spikenet);
% spikenet_ads.conditioning_type = 1;
% sim = spikenet_ads.start();  
% sim.run();  % This should produce output similar to
%             % `>> spikenet50_spike_trig()` 