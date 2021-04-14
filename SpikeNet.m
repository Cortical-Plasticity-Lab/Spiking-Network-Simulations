classdef SpikeNet
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
        folder
        subfolder = 'spikenet50'
        prefix = 'sn'
        mexname = 'spikenet50mex'
        conditioning_type
        stim_delay_ms
        stim_pulse_train
        stim_pulse_freq
        
    end
    
    methods
        function obj = untitled(inputArg1,inputArg2)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

