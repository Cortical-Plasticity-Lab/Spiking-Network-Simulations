classdef SimulationData < handle
   %SIMULATIONDATA A class for saving output data from simulations.
   
   properties
      activity          % activity time-series of each unit
      sweeps            % spike histogram counts for each unit
      sweep_count       % total number of sweeps
      LFP               % LFP separated by unit type and column
      evoked_potentials % Evoked potential (max of test range LFP - mean of baseline LFP)
      avetestLFP        % LFP averaging test result
      avetestLFP_sweeps % LFP averaging sweeps
      statest           % Spike-triggered averaging test results
      statest_sweeps    % Spike-triggered averaging sweeps
      psth              % Peri-stimulus time histogram
      psthn             % Number of sweeps in each histogram
      psthscale         % Time-scale of histogram
      psthedges         % Edges of histogram bins
      sta_sweeps
      scalems = [-50 100]
      nbins
      prebins
      sta
      scalex
      aveLFP
      aveLFP_sweeps
   end
   
   properties (Hidden,Access=public)
      itrial
      ntrials
      plotnames
      condnames
      LFPnames
      psthdefs
   end
   
   properties (Access=private)
      p                 % "Parent" SpikeNet
      foldername        % Name of output folder
   end
   
   methods
      function sim = SimulationData(p)
         %SIMULATIONDATA Class constructor
         %
         % Syntax:
         %  sim = SimulationData(p);
         %
         % Inputs:
         %  p - "Parent" SpikeNet handle.
         %
         % Output:
         %  sim - The SimulationData object.
         
         sim.p = p;
         sim.sweeps = zeros(p.time_steps, p.units); % spike histogram for each unit.
         sim.sweep_count = 0;
         sim.activity = zeros(p.units, p.time_steps, 'double'); % Save 10 seconds of activity for each unit.
         sim.LFP = zeros(p.LFP_count, p.time_steps, 'double');  % LFP separated by unit type and column.
         sim.evoked_potentials = zeros(9, 0);  % Evoked potential (max of test range LFP - mean of baseline LFP).
         
         sim.sta_sweeps = zeros(p.units,1);
         sim.nbins = (sim.scalems(2) - sim.scalems(1)) * p.msfs;
         sim.prebins = -sim.scalems(1) * p.msfs;
         sim.sta = zeros(p.units, sim.nbins);     % spike triggered averages
         sim.scalex = (-sim.prebins:sim.nbins-sim.prebins-1) / p.msfs;   % millisecond x-axis scale for averages
      end
      
      function run(sim)
         %RUN Run the sequence of trials
         
         % PSTH info
         npsth = length(sim.psthdefs);
         sim.psthscale = -100:0.5:100; % Scale in milliseconds
         sim.psthedges = sim.psthscale * sim.p.msfs; % Bin edges in samples
         sim.psth = zeros(length(sim.psthscale), npsth+10); % 10 extra PSTH for ST Ae1->Be and Ae1->Bi, or GT -> Ae/Ai/Be/Bi/Ce/Ci
         sim.psthn = zeros(1, npsth+10);                % Number of sweeps in each PSTH (used for converting counts to Hz)
         weightedges = (-sim.p.max_strength:10:sim.p.max_strength) + 5;
         
         % Sequence count info
         sequence_start_nocond = 1;
         sequence_start_test1 = sequence_start_nocond + sim.p.ntrials_off;
         sequence_start_cond = sequence_start_test1 + sim.p.ntrials_test;
         sequence_start_test2 = sequence_start_cond + sim.p.ntrials_on;
         sequence_end = sequence_start_test2 + sim.p.ntrials_test;
         sequence_length = sim.p.ntrials_off + sim.p.ntrials_test + ...
                           sim.p.ntrials_on + sim.p.ntrials_test;
         sim.ntrials = sim.p.ntrials;
         if sim.p.options.plot_PSTH || sim.p.options.plot_LFP
            main_fig = figure(...
               'Name', sprintf('Main Simulation-%02d Figure', sim.p.conditioning_type), ...
               'Position', [250 50 850 900], ...
               'Color', 'w');
         else
            main_fig = [];
         end
            
         for ii = 1:sim.ntrials
            sim.itrial = ii;
            sim.p.last_section_trial = 0;  % Flag to indicate last trial in the current sequence section.
            iseq = mod(sim.itrial - 1, sequence_length) + 1;
            if iseq == sequence_start_nocond
               sim.p.repeat_counter = sim.p.repeat_counter + 1;
               sim.p.conditioning_flag = 0;
               sim.p.train_on = 1;
               sim.p.sequence_counter = 1;
               section_label = 'Training no conditioning';
            elseif iseq == sequence_start_test1
               sim.p.conditioning_flag = 0;
               sim.p.train_on = 0;
               sim.p.sequence_counter = 1;
               section_label = 'Testing before conditioning';
            elseif iseq == sequence_start_cond
               sim.p.conditioning_flag = 1;
               sim.p.train_on = 1;
               sim.p.sequence_counter = 1;
               section_label = 'Training with conditioning';
               % Save test stim triggered LFP from previous section
               if sim.p.options.keep_sta
                  sim.avetestLFP = sim.aveLFP;
                  sim.avetestLFP_sweeps = sim.aveLFP_sweeps;
                  sim.statest = sim.sta;
                  sim.statest_sweeps = sim.sta_sweeps;
               end
            elseif iseq == sequence_start_test2
               sim.p.conditioning_flag = 0;
               sim.p.train_on = 0;
               sim.p.sequence_counter = 1;
               section_label = 'Testing after conditioning';
            else
               sim.p.sequence_counter = sim.p.sequence_counter + 1;
               nseq = iseq + 1; % Next trial sequence number
               if (nseq == sequence_start_test1) || (nseq == sequence_start_cond) || (nseq == sequence_start_test2) || (nseq == sequence_end)
                  sim.p.last_section_trial = 1;
               end
            end

            if sim.p.sequence_counter == 1  % Restart averages at the begining of each sequence section.
               sim.sta_sweeps = zeros(sim.p.units,1);
               sim.sta = zeros(sim.p.units, sim.nbins);
               sim.aveLFP = zeros(40,sim.nbins);
               sim.aveLFP_sweeps = 0;
               sim.psth = zeros(length(sim.psthscale), npsth+10);
               sim.psthn = zeros(1, npsth+10);
               sim.sweeps = zeros(sim.p.time_steps, sim.p.units);
               sim.sweep_count = 0;
               LFPstore = zeros(sim.p.n_cols, 0);
               sim.p.sumWeights_A_to_B = zeros(1,sim.p.time_steps); % Simulation sums weights A->B for specific training types.
            end

            if (sim.p.conditioning_type == 1)
               sim.p.rec_unit = sim.p.trigger_unit; % Spike triggered conditioning is on
            else
               sim.p.rec_unit = 0; % No spike trigger unit in other conditioning method
            end

            if (sim.p.conditioning_type == 4) || (sim.p.conditioning_type == 5) || (sim.p.conditioning_type == 6)
               % Tetanic conditioning is a fixed chance at each timestep
               chance = sim.p.tetanic_freq / sim.p.fs;
               tstamps = find(random('uniform', 0, 1, [1 floor(sim.p.fs * sim.p.conditioning_secs)]) < chance);
               short = find(diff(tstamps) < sim.p.stim_refractory);
               if ~isempty(short)
                  tstamps(short+1) = []; % Remove intervals that are shorter than our refractory period
               end
               sim.p.stim_source_times(:) = 0;
               sim.p.stim_source_times(tstamps) = 1;
               if sim.p.stim_pulse_train >= 2
                  % Pulse trains are at 300 Hz.
                  for ipulse = 1:sim.p.stim_pulse_train-1
                     sim.p.stim_source_times(tstamps + floor(3.3 * sim.p.msfs * ipulse)) = 1;
                  end
               end
            end

            % Run network for 1 iteration. Remember clock ticks of any conditioning
            % stimulations for later.
            sim.p.stim_output_times = zeros(sim.p.time_steps, 1, 'uint16'); % Will record the timesteps where conditioning stimuli were delivered.
            sim.p.trigger_times = zeros(sim.p.time_steps, 1, 'uint16');     % Clockticks where the conditioning triggers were detected (even if conditioning is turned off)
            [iUnit, tSpike] = spikenet50mex(sim.p.get_mex_struct(sim.itrial), sim.activity, sim.LFP);
            trig_timestamps = find(sim.p.trigger_times == 1);           % Convert to list of timestamps
            short = find(diff(trig_timestamps) < sim.p.stim_refractory);%   .. and remove short intervals (stim train intervals are 3.3 ms)
            if ~isempty(short)
               trig_timestamps(short+1) = [];
            end
            stim_timestamps = find(sim.p.stim_output_times == 1);
            clockticks = int32(stim_timestamps + (sim.itrial - 1) * sim.p.time_steps);
            sim.p.stim_output_clockticks = [sim.p.stim_output_clockticks; clockticks];
            for iu = 1:sim.p.units
               % Keep track of number of output spikes for each unit
               sim.p.unit_spike_counts(iu) = sim.p.unit_spike_counts(iu) + length(find(iUnit == iu));
            end

            % Calculate PSTHs
            % Trigger PSTHs are always calculated.  Group PSTHs are only calculated
            % for the final cycle (may need to extend this to several of the last
            % cycles for smoother graphs)
            for ipsth = 1:npsth
               psth_name = sim.psthdefs{ipsth}; % Text name of psth
               targ_name = psth_name(end-1:end);  % Target name (eg 'Ae')
               ref_name = psth_name(1:2);   % Reference name (eg 'Ae')

               % Get spike times for reference unit and target units
               if ref_name(2) == 'T'  % Special case for conditioning trigger times
                  ref_times = trig_timestamps;
               else
                  if (iseq >= sequence_start_test2) && (sim.p.ntrials_repeat == sim.p.repeat_counter)   % && (sim.p.last_section_trial == 1)
                     ref_times = tSpike(ismember(iUnit, find(strncmp(sim.p.unit_names, ref_name, 2))));
                  else
                     continue;
                  end
               end
               ref_times = ref_times(ref_times < 7 * sim.p.fs);
               targ_index = find(strncmp(sim.p.unit_names, targ_name, 2));
               targ_times = tSpike(ismember(iUnit, targ_index));

               % Sum a sweep for each reference spike time
               nref = length(ref_times);
               sim.psthn(ipsth) = sim.psthn(ipsth) + nref * length(targ_index);
               for iref = 1:nref
                  sweep_times = targ_times - ref_times(iref);
                  sweep_times = sweep_times((sweep_times >= sim.psthedges(1)) & (sweep_times <= sim.psthedges(end)));
                  if ~isempty(sweep_times)
                     sim.psth(:,ipsth) = sim.psth(:,ipsth) + histcounts(sweep_times, sim.psthedges);
                  end
               end
            end

            % Plot output sim.activity of representative units.

            for iun = 1:sim.p.units
               spikes = tSpike(iUnit == iun);
               if ~isempty(spikes)
                  sim.sweeps(spikes, iun) = sim.sweeps(spikes, iun) + 1; % Works because spikes should never overlap on the same index.
               end
            end

            if sim.p.options.plot_sweeps
               
               sim.sweep_count = sim.sweep_count + 1;
               for iplot = 1:6
                  iun = find(strncmp(sim.p.unit_names, sim.plotnames{iplot}, 2));
                  if ~isempty(iun)
                     subplot(6, 3, iplot * 3 - 2);
                     sumswp = sum(sim.sweeps(:, iun), 2) / length(iun);
                     sweep = 20 * sum(reshape(sumswp, round(sim.p.fs / 20), round(sim.p.time_steps / (sim.p.fs / 20)))) / sim.sweep_count;
                     plot(0:0.05:9.95, sweep);
                     xlabel([sim.plotnames{iplot} ' (sec)']);
                     ylabel('Hz');
                     ylim([0 50]);
                     if iplot == 1
                        if (sim.p.conditioning_type == 3)
                           titlestr = [sim.p.prefix ' ' datestamp ' (' sim.condnames{sim.p.conditioning_type+1} ' ' num2str(sim.p.stim_phase) ' deg, ' num2str(sim.p.stim_uV) ' uV)'];
                        elseif (sim.p.conditioning_type == 8)
                           band = [' ' num2str(sim.p.gamma_band(1)) '-'  num2str(sim.p.gamma_band(2)) ' Hz, '];
                           titlestr = [sim.p.prefix ' ' datestamp ' (' sim.condnames{sim.p.conditioning_type+1} band num2str(sim.p.stim_delay_ms) ' ms)'];
                        elseif (sim.p.conditioning_type == 0) || (sim.p.conditioning_type == 4) || (sim.p.conditioning_type == 5) || (sim.p.conditioning_type == 6)
                           titlestr = [sim.p.prefix ' ' datestamp ' (' sim.condnames{sim.p.conditioning_type+1} ')'];
                        else
                           titlestr = [sim.p.prefix ' ' datestamp ' (' sim.condnames{sim.p.conditioning_type+1} ' ' num2str(sim.p.stim_delay_ms) ' ms, ' num2str(sim.p.stim_uV) ' uV)'];
                        end
                        titlestr(titlestr == '_') = '-';
                        title(titlestr);
                     end
                  end
               end
            end

               
            % Plot spike triggered averages of LFP and a few other plots;
            if sim.p.options.plot_sta
               maxbin = sim.nbins - sim.prebins;

               % define which plots to show
               source(1) = find(strncmp(sim.p.unit_names, 'Ae', 2), 1); dest(1) = 1; % Ae1..n -> LFPA
               source(2) = source(1); dest(2) = 2; % LFPB
               source(3) = source(1); dest(3) = 3; % LFPC

               source(4) = find(strncmp(sim.p.unit_names, 'Be', 2), 1); dest(4) = 1; % Be1..n -> LFPA
               source(5) = source(4); dest(5) = 2; % LFPB
               source(6) = source(4); dest(6) = 3; % LFPC

               source(7) = find(strncmp(sim.p.unit_names, 'Ce', 2), 1); dest(7) = 1; % Ce1..n -> LFPA
               source(8) = source(7); dest(8) = 2; % LFPB
               source(9) = source(7); dest(9) = 3; % LFPC

               % plot spike triggered averages of LFP

               [b, a] = butter(1, [10 2500] / (sim.p.fs / 2)); %10Hz to 2500Hz butterworth filter
               for iplot = 1:9
                  iu = source(iplot);
                  iLFP = dest(iplot);

                  % Single unit spike triggered averages for first unit in each
                  % column.  Leave out seconds 8-10 which are involved in stimulus testing.
                  ts1 = tSpike(iUnit == iu);  % Spikes from first unit in source column (e.g Ae1 -> Col B)
                  LFPpos = filter(b, a, sim.LFP(iLFP, :));
                  LFPneg = filter(b, a, sim.LFP(iLFP+3, :));
                  for ispk = length(ts1):-1:1
                     t1 = ts1(ispk);
                     if (t1 > sim.nbins) && (t1 <= 7.5 * sim.p.fs - maxbin)
                        sim.sta_sweeps(iplot) = sim.sta_sweeps(iplot) + 1;
                        sim.sta(iplot, :) = sim.sta(iplot, :) + LFPpos((t1 - sim.prebins + 1):(t1 + maxbin));
                        sim.sta(iplot+50, :) = sim.sta(iplot+50, :) + LFPneg((t1 - sim.prebins + 1):(t1 + maxbin));
                     end
                  end

                  % Column wide spike triggered averages when training is off.
                  if (sim.p.train_on == 0)
                     ts1 = tSpike((iUnit >= iu) & (iUnit < iu + sim.p.n_excit)); % Spikes from all excitatory source units beginning with source name
                     for ispk = length(ts1):-1:1
                        t1 = ts1(ispk);
                        if (t1 > sim.nbins) && (t1 <= 75000 - maxbin)
                           sim.sta_sweeps(iplot+60) = sim.sta_sweeps(iplot+60) + 1;
                           sim.sta(iplot+60, :) = sim.sta(iplot+60, :) + LFPneg((t1 - sim.prebins + 1):(t1 + maxbin)) + LFPpos((t1 - sim.prebins + 1):(t1 + maxbin));
                        end
                     end
                  end

                  % Conditioning triggered averages for LFP Ae,Be,Ce,Ai,Bi,Ci, Aout,Bout,Cout
                  % column.  Leave out seconds 8-10 which are involved in stimulus testing.
                  ts1 = trig_timestamps((trig_timestamps > sim.nbins) & (trig_timestamps <= 7.5 * sim.p.fs - maxbin));  % Conditioning trigger events
                  if iplot <= 6
                     LFPpos = sim.LFP(iplot, :);
                  else
                     LFPpos = sim.LFP(iplot+3, :);
                  end
                  for ispk = length(ts1):-1:1
                     t1 = ts1(ispk);
                     sim.sta_sweeps(iplot+40) = sim.sta_sweeps(iplot+40) + 1;
                     sim.sta(iplot+40, :) = sim.sta(iplot+40, :) + LFPpos((t1 - sim.prebins + 1):(t1 + maxbin));
                  end

                  % Select correct plot
                  if iplot <= 6
                     subplot(6, 3, iplot * 3 - 1);
                  else
                     subplot(6, 3, (iplot - 6) * 3);
                  end

                  ave = sim.sta(iplot,:) / (1000 * sim.sta_sweeps(iplot)); % mV average
                  plot(sim.scalex, ave, 'Color', [0 0 0]);

                  hold on;
                  ave = sim.sta(iplot+50,:) / (1000 * sim.sta_sweeps(iplot));
                  plot(sim.scalex, ave, 'Color', [1 0 0]);

                  ave = (sim.sta(iplot,:) + sim.sta(iplot+50,:)) / (1000 * sim.sta_sweeps(iplot));
                  plot(sim.scalex, ave, 'Color', [0 0 1]);

                  ylim([-20 20]);
                  xlim([sim.scalems(1) sim.scalems(2)]);
                  name = sim.p.unit_names{iu};
                  xlabel([name ' -> LFP' sim.LFPnames{iLFP} ' sta (ms), n = ' num2str(sim.sta_sweeps(iplot))]);
                  ylabel('mV');
                  if iplot == 1
                     title('STA Unit Spike -> LFP');
                  end
                  if iplot == 7
                     title([section_label ' (' num2str(sim.itrial) ')']);
                  end
                  hold off;
               end

               % Update stimulus triggered LFP averages for test figures
               if (sim.p.train_on == 0)
                  sim.aveLFP_sweeps = sim.aveLFP_sweeps + 1;
                  r1 = sim.p.fs / 20; % Number of timesteps equal to 50 milliseconds
                  r2 = sim.p.fs / 40; % 25 milliseconds
                  iplot = 1;
                  if sim.p.last_section_trial
                     sim.evoked_potentials(:, end+1) = 0;
                  end
                  for icol = 1:3
                     LFPpos = filter(b, a, sim.LFP(icol, :));
                     LFPneg = filter(b, a, sim.LFP(icol+3, :));
                     LFPfilt = filter(b, a, sim.LFP(icol, :) + sim.LFP(icol+3, :));
                     for itime = 1:3
                        isample = 1 + 8 * sim.p.fs + (itime - 1) * 0.7 * sim.p.fs;  % Test stims are placed at specific times
                        range = round(isample - r1) : round(isample + sim.nbins - r1 - 1);
                        sim.aveLFP(iplot+20, :) = sim.aveLFP(iplot+20, :) + LFPpos(range); % Excite LFP
                        sim.aveLFP(iplot+30, :) = sim.aveLFP(iplot+30, :) + LFPneg(range); % Inhib LFP
                        sim.aveLFP(iplot, :) = sim.aveLFP(iplot, :) + LFPfilt(range); % Total LFP
                        sim.aveLFP(iplot+10, :) = sim.aveLFP(iplot+10, :) + sim.LFP(icol+7, range); % Synthetic EMG from output units
                        if sim.p.last_section_trial
                           % At the end of each testing section, calculate evoked
                           % potentials for each column to each other column
                           sim.evoked_potentials(iplot, end) = (max(sim.aveLFP(iplot,r1:r1+r2-1)) - mean(sim.aveLFP(iplot,r2:r1-1))) / sim.aveLFP_sweeps;
                        end
                        iplot = iplot + 1;
                     end
                  end

                  % Spike triggered averages for each unit in the Ae group onto LFPB
                  LFPfilt = filter(b, a, sim.LFP(2,:) + sim.LFP(2+3,:));
                  for iplot = 1:40
                     ts1 = tSpike(iUnit == iplot);  % Spikes from unit in source column (e.g Ae1 -> Col B)
                     for ispk = length(ts1):-1:1
                        t1 = ts1(ispk);
                        if (t1 > sim.nbins) && (t1 <= 7.5 * sim.p.fs - maxbin)
                           sim.sta_sweeps(iplot+70) = sim.sta_sweeps(iplot+70) + 1;
                           sim.sta(iplot+70, :) = sim.sta(iplot+70, :) + LFPfilt((t1 - sim.prebins + 1):(t1 + maxbin));
                        end
                     end
                  end
               end

            end

            % Calculate average weight from A -> B
            index = find((sim.p.weight_pre_unit <= sim.p.n_excit) & (sim.p.weight_post_unit > sim.p.n_excit + sim.p.n_inhib + sim.p.n_out) & (sim.p.weight_post_unit <= 2*sim.p.n_excit + 2*sim.p.n_inhib + sim.p.n_out));
            mean_weight_AB = mean(sim.p.weight_strength(index)) / sim.p.PSP_factor;
            index = find((sim.p.weight_pre_unit <= sim.p.n_excit) & (sim.p.weight_post_unit > 2*(sim.p.n_excit + sim.p.n_inhib + sim.p.n_out)));
            mean_weight_AC = mean(sim.p.weight_strength(index)) / sim.p.PSP_factor;

            % Display progress line
            progstr = ['Sweep ' num2str(sim.itrial) ' Spikes ' num2str(length(tSpike)) ' Stims ' num2str(sim.p.train_info(3)) ' WeightAB ' num2str(mean_weight_AB) ' WeightAC ' num2str(mean_weight_AC)];
            disp(progstr);
            fid = fopen([sim.p.output_folder '\progress.txt'], 'a');
            fprintf(fid, '%s\r\n', progstr);
            fclose(fid);

            if sim.p.conditioning_type == 1
               % For Spike triggered conditioning, plot PSTH of Ae1->Ae
               % Get spike times for reference unit and target units
               ref_times = tSpike(iUnit == sim.p.trigger_unit);
               %ref_times = tSpike(ismember(iUnit, find(strncmp(sim.p.unit_names, %'Ae1', 3), 1))); % Old version only supported Ae1 -> Be
               ref_times = ref_times(ref_times < sim.p.conditioning_secs * 7 * sim.p.fs); % Limit spikes to times before test stims
               targ_times = tSpike(ismember(iUnit, find(strncmp(sim.p.unit_names, 'Ae', 2))));
               ipsth = npsth+1; % The space we allocated for this plot in psth().

               % Sum a sweep for each reference spike time
               nref = length(ref_times);
               for iref = 1:nref
                  sweep_times = targ_times - ref_times(iref);
                  sweep_times = sweep_times((sweep_times >= sim.psthedges(1)) & (sweep_times <= sim.psthedges(end)));
                  if ~isempty(sweep_times)
                     sim.psth(:,ipsth) = sim.psth(:,ipsth) + histcounts(sweep_times, sim.psthedges);
                  end
               end
               
               if sim.p.options.plot_sta
                  subplot(6, 3, 12);
                  plot(sim.psthscale,  100 * sim.psth(:, ipsth) / sim.p.sequence_counter / sim.p.n_excit);
                  xlim([-50 50]);
                  trigname = sim.p.unit_names{sim.p.trigger_unit};
                  xlabel(['PSTH ' trigname ' -> Ae']);
                  ylabel('Hz');
               end

               % PSTH Ae1 -> all Be
               targ_times = tSpike(ismember(iUnit, find(strncmp(sim.p.unit_names, 'Be', 2))));
               ipsth = npsth+2; % The space we allocated for this plot in psth().

               % Sum a sweep for each reference spike time
               for iref = 1:nref
                  sweep_times = targ_times - ref_times(iref);
                  sweep_times = sweep_times((sweep_times >= sim.psthedges(1)) & (sweep_times <= sim.psthedges(end)));
                  if ~isempty(sweep_times)
                     sim.psth(:,ipsth) = sim.psth(:,ipsth) + histcounts(sweep_times, sim.psthedges);
                  end
               end
               
               if sim.p.options.plot_sta
                  subplot(6, 3, 15);
                  plot(sim.psthscale,  100 * sim.psth(:, ipsth) / sim.p.sequence_counter / sim.p.n_excit);
                  xlim([-50 50]);
                  xlabel(['PSTH ' trigname ' -> Be']);
                  ylabel('Hz');
               end

            elseif (sim.p.conditioning_type == 3) && sim.p.options.plot_LFP
               r2 = round(0.05 * sim.p.fs); % Graph 50 ms on either side of stim
               r1 = 2 * r2; % Timesteps equal to 100 ms
               sweep_ave = zeros(1,r1+1);
               sweep_index = find(sim.p.stim_output_times == 1);  % Find cycle triggers with full sweeps
               sweep_index = sweep_index((sweep_index > r2) & (sweep_index <  sim.p.time_steps - r2));
               n = length(sweep_index);
               if n > 0
                  % Beta band filtered LFPB is stored in LFP(7,:)
                  xaxis_ms = (-r2:r2) / sim.p.msfs;
                  for isweep = 1:n
                     sweep_start = sweep_index(isweep) - r2;
                     sweep_ave = sweep_ave + sim.LFP(7, sweep_start:sweep_start+r1);
                  end

                  subplot(6, 3, 12);
                  plot(xaxis_ms, sweep_ave / length(sweep_index));
                  xlabel(['Trig -> fLFPB (n = ' num2str(n) ')']);


                  sweep_index = find(sim.p.stim_target_times == 1); % find test triggers
                  sweep_ave = zeros(1,r1+1);
                  n = length(sweep_index);
                  if n > 0
                     for isweep = 1:n
                        sweep_start = sweep_index(isweep) - r2;
                        sweep_ave = sweep_ave + sim.LFP(7, sweep_start:sweep_start+r1);
                     end
                     subplot(6, 3, 15);
                     plot(xaxis_ms, sweep_ave / length(sweep_index));
                     xlabel(['Test -> fLFPB (n = ' num2str(n) ')']);

                  end
               end

            elseif ((sim.p.conditioning_type == 7) || (sim.p.conditioning_type == 8)) && sim.p.options.plot_LFP
               % EMG or Gamma triggered LFP-A average and ColA spikes
               r2 = round(0.05 * sim.p.fs); % Graph 50 ms on either side of triggers
               r1 = 2 * r2;  % Timesteps equal to 100 ms
               sweep_ave = zeros(1,r1+1);
               sweep_index = find(sim.p.stim_output_times == 1);  % Find gamma triggers with full sweeps
               sweep_index = sweep_index((sweep_index > r2) & (sweep_index <  sim.p.time_steps - r2));
               n = length(sweep_index);
               if n > 0
                  % Average LFPA (this is LFPA_pos + LFPA_neg)
                  xaxis_ms = (-r2:r2) / sim.p.msfs;
                  for isweep = 1:n
                     sweep_start = sweep_index(isweep) - r2;
                     sweep_ave = sweep_ave + sim.LFP(1, sweep_start:sweep_start+r1) + sim.LFP(4, sweep_start:sweep_start+r1);
                  end
                  subplot(6, 3, 12);
                  plot(xaxis_ms, sweep_ave / length(sweep_index));
                  xlabel(['T -> LFPA (n = ' num2str(n) ')']);

                  % EMG or Gamma filtered LFPA is returned in LFP(7,:)
                  sweep_ave = zeros(1,r1+1);
                  for isweep = 1:n
                     sweep_start = sweep_index(isweep) - r2;
                     sweep_ave = sweep_ave + sim.LFP(7, sweep_start:sweep_start+r1);
                  end
                  subplot(6, 3, 15);
                  plot(xaxis_ms, sweep_ave / length(sweep_index));
                  if sim.p.conditioning_type == 7
                     xlabel(['T -> EMGA (n = ' num2str(n) ')']);
                  else
                     xlabel(['T -> GammaA (n = ' num2str(n) ')']);
                  end
               end
            end

            if sim.p.options.plot_histograms
               % Replace last plot with a histogram of the weight distribution
               subplot(6, 3, 18);
               histogram(double(sim.p.weight_strength((sim.p.weight_strength ~= 0) & (sim.p.weight_training_rule ~= 0))) / sim.p.PSP_factor, weightedges);
               xlim([-sim.p.max_strength sim.p.max_strength]);
               ylim([0 500]);
               xlabel('Trainable weight distribution');
               drawnow;
            end
            
            % Accumulate LFP for a limited number of iterations
            % This is only for the PSD plot shown below the raster plot
            % Reduce the iterations allowed if space becomes an issue.
            if sim.p.options.plot_rasters && sim.p.sequence_counter <= 200
               for iLFP = 1:sim.p.n_cols
                  LFPcol = sim.LFP(iLFP,:) + sim.LFP(iLFP+3,:);
                  samples = length(LFPcol);
                  LFPstore(iLFP,end+1:end+samples) = LFPcol;
               end
            end

            % Save averages on last trial of each test section
            if (sim.p.last_section_trial == 1)
               meanWeightAB = sim.p.sumWeights_A_to_B ./ sim.p.PSP_factor ./ sim.p.sequence_counter;
               if (sim.p.train_on == 0)
                  if strcmp(section_label, 'Testing before conditioning')
                     if ~isempty(main_fig)
                        savefig(main_fig, [sim.foldername '\averages_pre_' num2str(sim.itrial) '.fig']);
                     end
                     save([sim.foldername '\meanWeightAB_pre.mat'], 'meanWeightAB');
                  else
                     if ~isempty(main_fig)
                        savefig(main_fig, [sim.foldername '\averages_post_' num2str(sim.itrial) '.fig']);
                     end
                     save([sim.foldername '\meanWeightAB_post.mat'], 'meanWeightAB');
                  end
               else
                  if sim.p.conditioning_flag
                     if ~isempty(main_fig)
                        savefig(main_fig, [sim.foldername '\averages_cond_' num2str(sim.itrial) '.fig']);
                     end
                     save([sim.foldername '\meanWeightAB_cond.mat'], 'meanWeightAB');
                  else
                     if ~isempty(main_fig)
                        savefig(main_fig, [sim.foldername '\averages_nocond_' num2str(sim.itrial) '.fig']);
                     end
                     save([sim.foldername '\meanWeightAB_nocond.mat'], 'meanWeightAB');
                  end
               end

               % Raster plot.  Plot spike times with LFPs and LFP Power spectrum.
               if sim.p.options.plot_rasters
                  raster_fig = figure('Name', sprintf('Simulation-%02d Rasters', sim.p.conditioning_type));
                  position = get(raster_fig, 'Position');
                  set(raster_fig,...
                     'Position', position .* [1 1 1.5 2], ...
                     'Color', 'w');
                  subplot(2,1,1);

                  samples = 1:floor(2*sim.p.fs); % Plot first 2 seconds
                  h = [0 0 0]; % Fill with handles to graphics appearing in the Legend.
                  colors = {'b', 'r', [0.7 0.7 0.5], 'c', 'm', 'y'};
                  LFPreserve = sim.p.n_excit;
                  ncolunits = sim.p.n_excit + sim.p.n_inhib + sim.p.n_out;
                  yreserve = LFPreserve + ncolunits + 20; % space to reserver for each column
                  ystart = 0;
                  yticks = [];
                  for iLFP = 1:sim.p.n_cols
                     colunits = (((iLFP - 1) * ncolunits) + 1) : (iLFP * ncolunits);
                     sindex = find((tSpike < samples(end)) & ismember(iUnit, colunits));
                     iu = iUnit(sindex);
                     miniu = min(iu);
                     iui = find(iu >= miniu + sim.p.n_excit + sim.p.n_inhib);
                     iu(iui) = iu(iui) + 20;
                     iui = find(iu > miniu + sim.p.n_excit);
                     iu(iui) = iu(iui) + 10;
                     iu = iu - colunits(1) + ystart + LFPreserve + 10;
                     plot(tSpike(sindex) ./ sim.p.fs, iu, 'k.');
                     hold on;
                     LFPcol = sim.LFP(iLFP,samples) + sim.LFP(iLFP+3,samples);
                     LFPcol = LFPcol - min(LFPcol);
                     LFPcol = (sim.p.n_excit * LFPcol / max(LFPcol)) + ystart;
                     h(iLFP) = plot((0:length(LFPcol)-1) ./ sim.p.fs, LFPcol, 'Color', colors{iLFP}, 'LineWidth', 2);
                     yticks(end+1) = ystart + floor(sim.p.n_excit/2);
                     yticks(end+1) = ystart + floor(sim.p.n_excit + sim.p.n_excit/2 + 10);
                     yticks(end+1) = ystart + floor(2*sim.p.n_excit + sim.p.n_inhib/2 + 20);
                     yticks(end+1) = ystart + floor(2*sim.p.n_excit + sim.p.n_inhib + sim.p.n_out/2 + 40);
                     ystart = ystart + yreserve + 40;
                  end
                  set(gca, 'fontsize', 16);
                  set(gca, 'fontweight', 'bold');
                  title(['Spikes: ' section_label ' (' num2str(sim.itrial) ')']);
                  xlabel('Time(Sec)');
                  set(gca, 'YTick', yticks);
                  set(gca, 'YTickLabel', {'LFPA', 'Ae', 'Ai', 'Am', 'LFPB', 'Be', 'Bi', 'Bm', 'LFPC', 'Ce', 'Ci', 'Cm'});
                  % legend(h, {'LFP A', 'LFP B', 'LFP C'});
                  xlabel('Time (Sec)');
                  ylim([0 max(iu)+10]);
                  xlim([0 samples(end)/sim.p.fs]); % Show first 1 second

                  subplot(2,1,2); %Spectral Density of LFP
                  freq = 0:1:199;
                  periodogram_type = 'Power'; % 'Power' or 'PSD'
                  for iLFP = 1:sim.p.n_cols
                     %[P, F] = periodogram(LFPstore(iLFP, :) - mean(LFPstore(iLFP, :)), [], freq, sim.p.fs, periodogram_type);
                     [P,F] = pwelch(LFPstore(iLFP, :) - mean(LFPstore(iLFP, :)),[],[],freq,sim.p.fs); % Usign pwelch over periodogram
                     plot(freq, P, 'Color', colors{iLFP}, 'LineWidth', 2);
                     %P = sum(reshape(P, 4, 50));
                     %plot((0:49)*4 + 2, P, 'Color', colors{iLFP}, 'LineWidth', 2);
                     hold on;
                  end
                  set(gca, 'fontsize', 16);
                  set(gca, 'fontweight', 'bold');
                  title('LFP Power Spectrum');
                  ylabel(periodogram_type);
                  xlabel('Hz');
                  legend({'LFP A', 'LFP B', 'LFP C'});
                  savefig(raster_fig, [sim.p.output_folder '\Raster-' num2str(sim.itrial) '.fig']);
                  close(raster_fig);
               end

               % Save PSTH figure for all T -> Column Groups
               if sim.p.options.plot_histograms
                  psth_fig = figure(...
                     'Name', 'PSTH Panels', ...
                     'Position', [865 50 850 900], ...
                     'Color', 'w' ...
                     );
                  for iplot = 1:9
                     % Plot PSTH
                     subplot(5,2,iplot);
                     data = sim.psth(:, iplot + 24) * 1000 / sim.psthn(iplot + 24) / (sim.psthscale(2) - sim.psthscale(1)); % Normalize to spikes per second
                     %bar(psthscale + 0.5 * (psthscale(2) - psthscale(1)), data, 1);  % Use this for solid fill bars
                     [xscale, yvals] = skyline(sim.psthscale, data);  % Create a skyline plot for the histogram
                     plot(xscale, yvals, 'LineWidth', 2);
                     ylabel('Hz');
                     ylim([0 60]);
                     xlim([-50 50]);
                     xticks(-50:10:50);
                     xticklabels({'-50', '', '', '', '', '0' , '', '', '', '', '50'});
                     xlabel([sim.psthdefs{iplot + 24} ', ' sim.LFPnames{iplot}]);
                     if iplot == 1
                        title(['PSTH ' char(titlestr)]);
                     end
                     if iplot == 2
                        title([section_label ' (' num2str(sim.itrial) ')']);
                     end

                     % Plot LFP of target column group
                     if iplot <= 6
                        ista = floor((iplot+1) / 2) + 40;  % Plot same LFP to both T->Ae and T->Ai graphs
                        data = (sim.sta(ista, :) + sim.sta(ista+3, :)) / 1000 / sim.sta_sweeps(ista); % Sum exitatory and inhibitory contributions to the LFP.
                     else
                        data = sim.sta(iplot+40, :) / 1000 / sim.sta_sweeps(iplot+40);
                     end
                     data = data - min(data);
                     hold on;
                     plot(sim.scalex, data, 'Color', [0.5,0.5,0.5], 'LineWidth', 2);

                  end
                  drawnow;
                  saveas(psth_fig, [sim.p.output_folder '\Cond_Trig_PSTH-' num2str(sim.itrial) '.fig']);
                  close(psth_fig);
               end
            end

         end
      end
      
      function fig = plot_PSTH(sim)
         %PLOT_PSTH Plot peri-stimulus time histograms
         fig = figure(...
            'Name', 'Peri-Stimulus Time Histogram', ...
            'Position', [10 50 850 920], ...
            'Color', 'w', ...
            'PaperOrientation', 'Portrait', ...
            'PaperType', 'usletter');
         iplot = 1;
         for icol = 1:3
             for ipsth = 1:8
                 subplot(8,3, (ipsth -1) * 3 + icol);
                 vals = sim.psth(:, iplot) * 1000 / sim.psthn(iplot) / (sim.psthscale(2) - sim.psthscale(1));
                 % zerobin = find(sim.psthedges == 0);
                 % vals(zerobin) = 0.5 * (vals(zerobin-1) + vals(zerobin+1)); % average over center bin.
                 plot(sim.psthscale, vals);
                 xlim([-25 25]);
                 ylim([0 25]);
                 xlabel(sim.psthdefs{iplot});
                 if iplot == 1
                     titlestr2 = ['/' sim.p.subfolder '/' sim.p.prefix ' PSTH ' datestamp];
                     titlestr2(titlestr2 == '_') = '-';
                     title(titlestr2);
                 end
                 iplot = iplot + 1;
                 hold on;
             end
         end
         drawnow
         saveas(fig, [sim.foldername '\PSTH.fig']);
      end
   end
   
   methods (Access=private)
      function init_defaults_(sim)
         sim.foldername = sim.p.output_folder;
         sim.plotnames = {'Ae'; 'Be'; 'Ce'; 'Ao'; 'Bo'; 'Co'};
         sim.LFPnames = {'A'; 'B'; 'C'; 'A-'; 'B-'; 'C-'; 'A0'; 'B0'; 'C0'};
         sim.condnames = {'No Conditioning'; 'ST'; 'PP'; 'CT'; 'Tetanic-A'; 'Tetanic-B'; 'Tetanic-C'; 'ET'; 'GT'; 'AT'};
         sim.psthdefs = { ...
             'Ae -> Ae'; 'Ae -> Ai'; 'Ae -> Be'; 'Ae -> Bi'; 'Ae -> Ce'; 'Ae -> Ci'; 'Ai -> Ae'; 'Ai -> Ai'; ...
             'Be -> Ae'; 'Be -> Ai'; 'Be -> Be'; 'Be -> Bi'; 'Be -> Ce'; 'Be -> Ci'; 'Bi -> Be'; 'Bi -> Bi'; ...
             'Ce -> Ae'; 'Ce -> Ai'; 'Ce -> Be'; 'Ce -> Bi'; 'Ce -> Ce'; 'Ce -> Ci'; 'Ci -> Ce'; 'Ci -> Ci'; ...
             ' T -> Ae'; ' T -> Ai'; ' T -> Be'; ' T -> Bi'; ' T -> Ce'; ' T -> Ci'; ' T -> Ao'; ' T -> Bo'; ' T -> Co'}; 
      end
   end
   
end
   