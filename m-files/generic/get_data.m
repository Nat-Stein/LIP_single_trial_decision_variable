function [H, unitIdxLIP_Tin, unitIdxLIP_Tout, t, RT, include, coh, direction, choice, correct, ...
    unitIdxLIP_dotsInRF, pulse_on, pulse_off, pulse_size] = get_data(dt_ms, neurons_to_include, dataset,par,varargin)

% AZ 2022 modified by NST 2023

dt_nans_before_RT = 0.05; % in seconds

for i=1:2:length(varargin)
    if isequal(varargin{i},'dt_rel_RT')
        dt_nans_before_RT = varargin{i+1};
    end
end

extract_neural_data_flag = 1;
if isnan(dt_ms)
    extract_neural_data_flag = 0;
end

if nargin==0 || isempty(dt_ms)
    dt_ms = 5;
end

if nargin<3 || isempty(dataset)
    dataset = 1;
end

% --------
% NSt edit:
datadir = par.rawPath;
filename = par.recordingName;
load(fullfile(datadir,filename));

unitIdxLIP_Tin = par.unit_Tin;
unitIdxLIP_Tout = par.unit_Tout;
unitIdx_DinRFcC = par.unit_MinC;
unitIdx_DinRFiC = par.unit_MinI;
unitIdxLIP_dotsInRF = par.unit_MinC;
minCells = struct('DinRFc', unitIdx_DinRFcC ,'DinRFi',unitIdx_DinRFiC);
% --------

%%

if nargin<2 || isempty(neurons_to_include)
    neurons_to_include = 1:length(d(1).spCellPop);
end

%%

trialType = [d.trialType]';

idx_dots = trialType==20; % dots task

dotsOn = [d.dotsOn]';
dotsOff = [d.dotsOff]';
saccadeDetected = [d.saccadeDetected]';

RT = saccadeDetected - dotsOn;
% include = idx_dots == 1 & ~isnan(dotsOn) & ~isnan(dotsOff) & ~isnan(RT);
complete = [d.complete]'; % NSt
include = idx_dots == 1 & ~isnan(dotsOn) & ~isnan(dotsOff) & ~isnan(saccadeDetected) & complete;% NSt

targetsOn = [d.targetsOn]';

coh = [d.sCoh]';
direction = [d.sDir]';
choice = [d.choice]';
correct = [d.correct]';


if extract_neural_data_flag
    nTr = length(dotsOn);
    pre_t = 0.2;
    post_t = 5; % 5 seconds max
    nTimeSteps = round(1000*(pre_t+post_t)/dt_ms) + 1;
    nNeurons = length(neurons_to_include);
    H = nan(nNeurons, nTimeSteps, nTr);
    R_ton = nan(nNeurons, nTr);
    for iTr = 1:nTr
        if ~isnan(dotsOn(iTr))
            tini = dotsOn(iTr) - pre_t;
            tend = dotsOn(iTr) + post_t;
            spk = d(iTr).spCellPop(neurons_to_include);
            dt = dt_ms/1000;
            [t,H(:,:,iTr)] = spikeanalysis.spk_to_hist(spk,tini,tend,dt);
            t = t - dotsOn(iTr);
            H(:,(t+dt/2)>[saccadeDetected(iTr) - dotsOn(iTr) - dt_nans_before_RT],iTr) = nan; % NaN spikes after RT

            % target-window
            wt = 0.1; % sec
            [~,aux] = spikeanalysis.spk_to_hist(spk,targetsOn(iTr) - wt, targetsOn(iTr), wt);
            R_ton(:,iTr) = aux(:,1);

        end
    end

    % get the spike times - new
    spike_times_ms = cell(nTr, nNeurons);
    for iTr = 1:nTr
        spk = d(iTr).spCellPop(neurons_to_include);
        for j=1:nNeurons
            spike_times_ms{iTr,j} = 1000*(spk{j}-dotsOn(iTr));
        end
    end



    % extract a higher-res version (1ms), with just the Tin neurons
    dt_high = 1/1000;
    nTimeStepsHigh = [post_t + pre_t]/dt_high + 1;
    DV = nan(nTr, nTimeStepsHigh);
    for iTr = 1:nTr
        if ~isnan(dotsOn(iTr))
            tini = - pre_t;
            tend = + post_t;
            spk = d(iTr).spCellPop(unitIdxLIP_Tin);
            % group them
            spk = {cat(1,spk{:}) - dotsOn(iTr)};
            [td,DV(iTr,:)] = spikeanalysis.spk_to_hist(spk,tini,tend,dt_high);
            %DV(iTr,(td+dt/2)>[saccadeDetected(iTr) - dotsOn(iTr) - 0.05],iTr) = nan; % NaN spikes after RT

        end
    end
    % remove post-decision samples
    DV = motionenergy.remove_post_decision_samples(DV,td, RT); % set nan's after RT
    
    % divide by the number of neurons grouped
    DV = DV / length(unitIdxLIP_Tin);

    DV = DV / (td(2)-td(1)); % for units of firing rate

    %smoothed high
    ns = round(0.08/(td(2)-td(1))); % smoothing time steps for 80ms
    dv_tin_high_smooth = conv2(DV, ones([1,ns]),'same')/ns;
    dv_tin_high_smooth(:,td<-0.1) = nan;
    dv_tin_high_smooth = motionenergy.remove_post_decision_samples(dv_tin_high_smooth, td, RT-dt_nans_before_RT);
    
    
else
    t = [];
    H = [];
    R_ton = [];
    
end

pulse_on = [d.pulseOn]' - [d.dotsOn]';
pulse_off = [d.pulseOff]' - [d.dotsOn]';
pulse_size = [d.pulseSize]';



if nargout==1
   I = include & ~isnan(RT) & [d.complete]'==1;
   if extract_neural_data_flag
       out = struct('H',H(:,:,I), 'unitIdxLIP_Tin',unitIdxLIP_Tin, 'unitIdxLIP_Tout', unitIdxLIP_Tout, ...
           't',t, 'RT',RT(I), 'included',I, 'coh',coh(I), 'direction',direction(I), 'choice',choice(I), ...
           'correct',correct(I), 'unitIdxLIP_dotsInRF',unitIdxLIP_dotsInRF, 'pulse_on',pulse_on(I), ...
           'pulse_off',pulse_off(I), 'pulse_size',pulse_size(I),'dataset',filename,...
           'R_ton',R_ton(:,I),'td',td,'DV',DV(I,:),'DV_smooth',dv_tin_high_smooth(I,:),...
           'nneurons',size(H,1), 'ntimes',size(H,2),'ntrials',sum(I));
       
       out.spike_times_ms = spike_times_ms(I,:);

   else
       out = struct('RT',RT(I), 'included',I, 'coh',coh(I), 'drection',direction(I), 'choice',choice(I), ...
           'correct',correct(I), 'pulse_on',pulse_on(I), ...
           'pulse_off',pulse_off(I), 'pulse_size',pulse_size(I),'dataset',filename);
   end
   
   
   file_new_tin = ['../prepro_recalc_Tin_Tout/d',num2str(dataset),'.mat'];
   if extract_neural_data_flag && exist(file_new_tin,'file')
       aux = load(file_new_tin);
       out.idx_Tin = aux.idx_Tin;
       out.idx_Tout = aux.idx_Tout;
   end
   
   out.minCells = minCells;
   
   H = out; % first output

end


end

