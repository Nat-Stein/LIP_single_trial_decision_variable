function [spk, spkMap, spkView, par] = reformDat(par, dd)
%% Generate spiking histograms for all trials and neurons
% Inputs
%   - par = basic information about the recording session and parameters
%           including 
%           the time vectors par.tsl (time aligned to motion onset in ms)
%           and par.trl (time aligned to motion onset in ms)
%   - dd  = data file of the recording session
% Outputs
%   - spk = matrix of size [neurons x trials x time points] specifying
%           whether neuron n spiked (1) or not (0) during trial tt at time 
%           point t during the RDM discrimination task
%           spk.SL = aligned to motion onset
%           spk.SLnan = aligned to motion onset with NaN starting 100ms
%           before saccade on every trial
%           spk.RL = aligned to saccade
%           spk.RLnan = aligned to saccade with NaN for time points before 200ms
%           after motion onset on every trial
%           spk.TRG = aligned to reponse target onset
%           spk.TRGnan = aligned to response target onset with NaN starting 
%           at the time the motion stimulus apears on every trial
%  - spkMap = same as spk, but during the instructed saccade task (mapping
%           task)
%  - spkView = same as spk, but during the passive motion viewing task


beh = par.beh;
numTrD = length(beh.datMat2Cell);
behMap = par.behMap;
numTrM = length(behMap.datMat2CellMap);
behView = par.behView;
numTrV = length(behView.dot_on);

numCells = length(dd(1).spCellPop);
par.numCells = numCells;

spk.SL = nan(numCells, numTrD, length(par.tsl));
spk.RL = nan(numCells, numTrD, length(par.trl));
spk.SLnan = nan(numCells, numTrD, length(par.tsl));
spk.RLnan = nan(numCells, numTrD, length(par.trl));
spk.TRG = nan(numCells, numTrD, length(par.trg));
spk.TRGnan = nan(numCells, numTrD, length(par.trg));

% Set-up spiking matrices
spkMap.mat = nan(numCells, numTrM,5000);
spkMap.SL = nan(numCells, numTrM, length(par.tsl));
spkMap.RL = nan(numCells, numTrM, length(par.trl));
spkMap.FP = nan(numCells, numTrM, length(par.tfp));
spkMap.SLnan = nan(numCells, numTrM, length(par.tsl));
spkMap.RLnan = nan(numCells, numTrM, length(par.trl));
spkMap.FPnan = nan(numCells, numTrM, length(par.tfp));

spkView.SL = nan(numCells, numTrV, length(par.tsl));
spkView.SLnan = nan(numCells, numTrV, length(par.tsl));

% ---------------------------------------------
% RDM discrimination task start

for dc = 1 : numCells
    
    % disp(num2str(dc))
    % Dot trials
    if numTrD > 0
        
        
        for t2 = 1 : numTrD
            
            tt = beh.datMat2Cell(t2);
            
            t_doton = beh.dot_on(t2);
            t_dotoff = beh.dot_off(t2);
            t_targPres = beh.t_targPres(t2);
            t_respMad = beh.saccade(t2);
            
            these_spikes = dd(tt).spCellPop{dc};
            
            isisTRL(t2, :) = hist(diff(these_spikes), par.isihistbins);
            numSpikesRTL(t2) = length(these_spikes);
            ISIviolationsTRL(t2) = length(find(diff(these_spikes) < par.minISI));
            
            % Spike trains aligned to target onset (2 targets)
            histbins = [t_targPres + par.eplims{3}(1) : par.dt : t_targPres + par.eplims{3}(2)];
            thisTRL = these_spikes(find(these_spikes >= histbins(1) & these_spikes <= histbins(end)));
            spk.TRG(dc, t2,1:length(histbins)) = hist(thisTRL, histbins);
            
            % Spike trains aligned to target onset (2 targets) - nan after response
            thisDotDel = t_doton - t_targPres;
            spk.TRGnan(dc, t2, :) = spk.TRG(dc, t2, :);
            spk.TRGnan(dc, t2, find(histbins > t_doton)) = NaN;
            
            % Stimulus-locked spike trains
            histbins = [t_doton + par.eplims{1}(1) : par.dt : t_doton + par.eplims{1}(2)];
            thisTRL = these_spikes(find(these_spikes >= histbins(1) & these_spikes <= histbins(end)));
            spk.SL(dc, t2,1:length(histbins)) = hist(thisTRL, histbins);
            
            % Stimulus-locked spike trains - nan after response
            thisRT = beh.rt(t2);
            spk.SLnan(dc, t2, :) = spk.SL(dc, t2, :);
            spk.SLnan(dc, t2, find(par.tsl > thisRT - par.cutSL)) = NaN;
            
            % Resonse-locked spike trains
            histbins = [t_respMad + par.eplims{2}(1) : par.dt : t_respMad + par.eplims{2}(2)];
            thisTRL = these_spikes(find(these_spikes >= histbins(1) & these_spikes <= histbins(end)));
            spk.RL(dc, t2,1:length(histbins)) = hist(thisTRL, histbins);
            
            % Response-locked spike trains - nan before dot onset
            spk.RLnan(dc, t2, :) = spk.RL(dc, t2, :);
            spk.RLnan(dc, t2, find(par.trl < - beh.rt(t2) + par.cutRL)) = NaN;
            
        end
        spk.isis{dc} = nansum(isisTRL, 1);
        spk.numSpikes(dc) = nansum(numSpikesRTL);
        spk.ISIviolations(dc) = nansum(ISIviolationsTRL);
    end
    
    % RDM discrimination task end
    % ---------------------------------------------
    % Mapping trials start
    
    if numTrM > 0
        
        % -------------------------------------------------
        
        clear isisTRL numSpikesRTL ISIviolationsTRL
        
        for t3 = 1 : numTrM
            
            tt = behMap.datMat2CellMap(t3); % Not all of these have all of the following events!
            
            t_fpOn = behMap.t_fpOn(t3);
            t_targPres = behMap.t_targPres(t3);
            t_fpOff = behMap.t_fpOff(t3);
            t_response = behMap.t_response(t3);
            
            behMap.targDelay(t3) = t_targPres - t_fpOn;
            
            these_spikes = dd(tt).spCellPop{dc};
            
            isisTRL(t3, :) = hist(diff(these_spikes), par.isihistbins);
            numSpikesRTL(t3) = length(these_spikes);
            ISIviolationsTRL(t3) = length(find(diff(these_spikes) < par.minISI));
            
            % Spike trains aligned to fixation onset
            histbins = [t_fpOn + par.eplims{3}(1) : par.dt : t_fpOn + par.eplims{3}(2)];
            thisTRL = these_spikes(find(these_spikes >= histbins(1) & these_spikes <= histbins(end)));
            spkMap.FP(dc, t3,1:length(histbins)) = hist(thisTRL, histbins);
            
            % Spike trains aligned to fixation onset - nan after target onset
            thisCut = behMap.targDelay(t3);
            spkMap.FPnan(dc, t3, :) = spkMap.FP(dc, t3, :);
            spkMap.FPnan(dc, t3, find(par.tfp > thisCut - par.cutFPoff)) = NaN; % is / 1000 for dot trials
            
            % ----------------------------------------------
            % Target-aligned spike trains
            histbins = [t_targPres + par.eplims{1}(1) : par.dt : t_targPres + par.eplims{1}(2)];
            thisTRL = these_spikes(find(these_spikes >= histbins(1) & these_spikes <= histbins(end)));
            spkMap.SL(dc, t3,1:length(histbins)) = hist(thisTRL, histbins);
            
            % Stimulus-locked spike trains - nan after fixation off
            thisCut = behMap.fixHold(t3);
            spkMap.SLnan(dc, t3, :) = spkMap.SL(dc, t3, :);
            spkMap.SLnan(dc, t3, find(par.tsl > thisCut - par.cutFPoff)) = NaN;  % is / 1000 for dot trials
            
            % Resonse-locked spike trains
            histbins = [t_response + par.eplims{2}(1) : par.dt : t_response + par.eplims{2}(2)];
            thisTRL = these_spikes(find(these_spikes >= histbins(1) & these_spikes <= histbins(end)));
            spkMap.RL(dc, t3,1:length(histbins)) = hist(thisTRL, histbins);
            
            % Response-locked spike trains - nan before target onset
            thisCut = behMap.fixHold(t3);
            spkMap.RLnan(dc, t3, :) = spkMap.RL(dc, t3, :);
            spkMap.RLnan(dc, t3, find(par.trl < - thisCut + par.cutRLMap)) = NaN;
        end
        spkMap.isis{dc} = nansum(isisTRL, 1);
        spkMap.numSpikes(dc) = nansum(numSpikesRTL);
        spkMap.ISIviolations(dc) = nansum(ISIviolationsTRL);
    end
    
    % Mapping trials end
    % ------------------------------------------------------------------
    % Motion viewing trials start
    
    if numTrV > 0
        
        clear isisTRL numSpikesRTL ISIviolationsTRL
        
        for t4 = 1 : numTrV
            
            tt = behView.datMat4Cell(t4);
            
            t_doton = behView.dot_on(t4);
            t_dotoff = behView.dot_off(t4);
            t_targPres = behView.t_targPres(t4);
            t_respMad = behView.saccade(t4);
            
            these_spikes = dd(tt).spCellPop{dc};
            
            isisTRL(t4, :) = hist(diff(these_spikes), par.isihistbins);
            numSpikesRTL(t4) = length(these_spikes);
            ISIviolationsTRL(t4) = length(find(diff(these_spikes) < par.minISI));
            
            % Stimulus-locked spike trains
            histbins = [t_doton + par.eplims{1}(1) : par.dt : t_doton + par.eplims{1}(2)];
            thisTRL = these_spikes(find(these_spikes >= histbins(1) & these_spikes <= histbins(end)));
            spkView.SL(dc, t4,1:length(histbins)) = hist(thisTRL, histbins);
            
            % Stimulus-locked spike trains - nan after response
            thisRT = behView.dot_off(t4) - behView.dot_on(t4);
            spkView.SLnan(dc, t4, :) = spkView.SL(dc, t4, :);
            spkView.SLnan(dc, t4, find(par.tsl > thisRT - par.cutSL)) = NaN;
        end
        spkView.isis{dc} = nansum(isisTRL, 1);
        spkView.numSpikes(dc) = nansum(numSpikesRTL);
        spkView.ISIviolations(dc) = nansum(ISIviolationsTRL);
    end
end


end

