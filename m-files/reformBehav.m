function [beh, behMap, behView] = reformBehav(par, dd)
% reformBehav reformat behavioral indicators per trial into structures for
% different trial types
% Inputs
%       par = structure defining file locations and a few additional
%           paprameters
%       dd = structure including behavioral indicators per trial
% Outputs
%       beh = structure with behavioral indicators for random dot motion
%               discrimination trials
%       behMap = structure with behavioral indicators for mapping trials
%       behView = structure with behavioral indicators for random dot motion
%               viewing trials


% Extract triggers and behavioural indicators
parTrigList

% Define task indicators
dotIndSacc = 20;     % motion discrimination with saccadic response
dotIndReach = 21;    % motion discrimination with reach response
orInd = 1;           % overalp reaches
mrInd = 11;          % memory reaches
osInd2 = 22;         % overlap saccades

if par.date(end) ~= 'N'
    osInd = 2;      % overlap saccades
    msInd = 3;      % memory saccades
    viewInd = 50;   % motion viewing
else
    osInd = 0;      % overlap saccades
    msInd = 10;      % memory saccades
    viewInd = 24;   % motion viewing
end

for tt = 1 : length(dd)
    trialType(tt) = dd(tt).trialType;
end


% Sort behavioural indicators into out.dataMat (dots) and out.dataMap (mapping) matrices
t2 = 0; t3 = 0;
beh = [];
behMap = [];
for tt = 1 : length(dd)
    
    % ============================================================
    % MOTION DISCRIMINATION TASK
    if dd(tt).trialType == dotIndSacc %|| dd(tt).trialType == dotIndReach
        
        if isnan(dd(tt).fixBreakTime) && dd(tt).complete == 1
            t2 = t2 + 1;
            
            beh.datMat2Cell(t2) = tt;

            beh.dot_on(t2) = dd(tt).dotsOn;             % Time stamp of motion onset
            beh.dot_off(t2) = dd(tt).dotsOff;           % Time stamp of motion offset
            beh.saccade(t2) = dd(tt).saccadeDetected;   % Time stamp of saccade onsets
            beh.dot_dir(t2) = dd(tt).direction;         % Direction of dot motion in degrees
            beh.dot_coh(t2) = dd(tt).coh;               % Absolute coherence
            beh.sig_coh(t2) = dd(tt).sCoh;              % Signed coherence
            beh.seed_base(t2) = dd(tt).seedBase;        % Seed base of random number generator
            beh.seed_var(t2) = dd(tt).seedVar;          % Seed var of random number generator
            beh.correct(t2) = dd(tt).correct;           % 1 = correct response, 0 = error
            beh.cho_trg(t2) = dd(tt).choice;            % 0 = left, 1 = right
            beh.rt(t2) = dd(tt).saccadeDetected - dd(tt).dotsOn;    % Reaction time in s
            beh.dot_dur(t2) = dd(tt).dotsOff - dd(tt).dotsOn;       % Duration dots were displayed on screen in s
            beh.dot_sp(t2) = dd(tt).dotsSpeed;          % Speed of random dot motion in degrees of visual angle/s * 10
            beh.dot_apr(t2) = dd(tt).dotsDiameter;      % Diameter of random dot motion patch in degrees of visual angle * 10
            beh.task(t2) = dd(tt).trialType - 19;       % 1 = saccade task, 2 = reach task
            beh.effector(t2) = 1;                       % Always saccade response
            beh.targ_acq(t2) = dd(tt).saccadeComplete;  % Time stamp when gaze landed on response target
            beh.fix_x(t2) = dd(tt).fixPos(1);           % x position of fixation point
            beh.fix_y(t2) = dd(tt).fixPos(2);           % y position of fixation point
            beh.targ1_x(t2) = dd(tt).targ1Pos(1);       % x position of target 1
            beh.targ1_y(t2) = dd(tt).targ1Pos(2);       % y position of target 1
            beh.targ2_x(t2) = dd(tt).targ2Pos(1);       % x position of target 2
            beh.targ2_y(t2) = dd(tt).targ2Pos(2);       % y position of target 2
            beh.t_targPres(t2) = dd(tt).targetsOn;      % Time stamp when response targets were presented
            beh.pulseTrial(t2) = dd(tt).pulseTrial;     % 1 = pulse was presented on this trial, 0 = no pulse
            
            beh.t_pulseon(t2) = dd(tt).pulseOn;         % Time
            beh.t_pulseoff(t2) = dd(tt).pulseOff;
            beh.pulseSize(t2) = dd(tt).pulseSize;
            beh.pulseDur(t2) = dd(tt).pulseDuration;
            if isnan(dd(tt).pulseSize)
                beh.pulseCoh(t2) = dd(tt).sCoh;
            else
                beh.pulseCoh(t2) = dd(tt).pulseSize/1000 + dd(tt).sCoh;
            end
            beh.trialid(t2) = dd(tt).trialType;
            beh.rewDur(t2) = dd(tt).rewardSize;
            % choiceG(t2) = dd(tt).choice;
        end
        
    % ============================================================
    % MAPPING TASK
    elseif dd(tt).trialType == osInd || dd(tt).trialType == msInd ||...
            dd(tt).trialType == orInd || dd(tt).trialType == mrInd % Mapping trial (add osInd2?)
        if isnan(dd(tt).fixBreakTime)
            t3 = t3 + 1;
            behMap.datMat2CellMap(t3) = tt;
            
            behMap.t_fpOn(t3) = dd(tt).fixOn;
            behMap.t_targPres(t3) = dd(tt).targetsOn;
            
            behMap.t_fpOff(t3) = dd(tt).fixOff;
            behMap.t_response(t3) = dd(tt).saccadeDetected; 
            behMap.t_targ_acq(t3) = dd(tt).saccadeComplete; 
            if (dd(tt).trialType == orInd || dd(tt).trialType == mrInd) &&...
                    isnan(dd(tt).saccadeDetected) && ~isnan(dd(tt).saccadeComplete)
                behMap.t_response(t3) = behMap.t_fpOff(t3)+0.25; 
            end
            behMap.targDelay(t3) = dd(tt).targetsOn - dd(tt).fixOn;
            if dd(tt).trialType == osInd
                behMap.targDur(t3) = dd(tt).rewardOn - dd(tt).targetsOn;
                behMap.t_targOff(t3) = dd(tt).rewardOn;
            else
                behMap.targDur(t3) = 0.2;          % duration memory saccade target is displayed
                behMap.t_targOff(t3) = dd(tt).targetsOn + 0.2;
            end
            behMap.fixHold(t3) = behMap.t_response(t3) - dd(tt).fixAcquired; 
            behMap.goRT(t3) = behMap.t_response(t3) - dd(tt).fixOff; 
            behMap.respDur(t3) = dd(tt).saccadeComplete - behMap.t_response(t3); 
            behMap.targ_x(t3) = dd(tt).targ1Pos(1) / 10;
            behMap.targ_y(t3) = dd(tt).targ1Pos(2) / 10;
            behMap.task(t3) = rem(dd(tt).trialType, 10) + 1;
            behMap.effector(t3) = behMap.task(t3);
            if par.date(end) ~= 'N'
                behMap.overMem(t3) = dd(tt).trialType - 1;
            else
                behMap.overMem(t3) = dd(tt).trialType/10 + 1;
            end
            behMap.fix_x(t3) = dd(tt).fixPos(1) / 10;
            behMap.fix_y(t3) = dd(tt).fixPos(2) / 10;
            
            behMap.trialid(t3) = dd(tt).trialType;
            behMap.trg_ct(t3) = 1; % target count: number of target within trial
            behMap.trlNum(t3) = t3; % trial number (ignoring target count)
            behMap.t_fpOnEye(t3) = NaN;
            behMap.t_targPresEye(t3) = NaN;
            behMap.t_targOffEye(t3) = NaN;
            behMap.t_fpOffEye(t3) = NaN;
            behMap.t_responseEye(t3) = NaN;
            behMap.t_targ_acqEye(t3) = NaN;
        end
    end
end
numTr = length(beh.datMat2Cell);
beh.dataCell = cell(numTr, 2);
beh.stimCell = cell(numTr, 2);
beh.reachPath = cell(numTr, 2);


% ============================================================
% MOTION VIEWING TASK

if par.date(end) ~= 'N'
    viewInd = 50;   % motion viewing
else
    viewInd = 24;   % motion viewing
end

behView = []; t4 = 0; 
for tt = 1 : length(dd)
    
    if dd(tt).trialType == viewInd % motion viewing
        
        validTrial = 0; 
        if par.date(end) ~= 'N'
            validTrial = isnan(dd(tt).fixBreakTime) & dd(tt).complete == 1 & ~isnan(dd(tt).dotsOff);
        else
            validTrial = isnan(dd(tt).fixBreakTime) & ~isnan(dd(tt).dotsOff) & isnan(dd(tt).saccadeDetected);
        end
        
        if validTrial
            t4 = t4 + 1;
            
            behView.datMat4Cell(t4) = tt;
            
            behView.dot_on(t4) = dd(tt).dotsOn;
            behView.dot_off(t4) = dd(tt).dotsOff;
            behView.saccade(t4) = dd(tt).saccadeDetected;
            behView.dot_dir(t4) = dd(tt).direction;
            behView.dot_coh(t4) = dd(tt).coh;
            behView.sig_coh(t4) = dd(tt).sCoh;
            behView.seed_base(t4) = dd(tt).seedBase;
            behView.seed_var(t4) = dd(tt).seedVar;
            behView.correct(t4) = dd(tt).correct;
            behView.cho_trg(t4) = dd(tt).choice;
            behView.rt(t4) = (dd(tt).saccadeDetected - dd(tt).dotsOn) * 1000;
            behView.dot_dur(t4) = dd(tt).dotsOff - dd(tt).dotsOn;
            behView.dot_sp(t4) = dd(tt).dotsSpeed;
            behView.dot_apr(t4) = dd(tt).dotsDiameter;
            behView.task(t4) = 1;                      % Always saccade task
            behView.effector(t4) = 1;                  % Response always with saccade
            behView.targ_acq(t4) = dd(tt).saccadeComplete;
            behView.fix_x(t4) = dd(tt).fixPos(1);
            behView.fix_y(t4) = dd(tt).fixPos(2);
            behView.targ1_x(t4) = dd(tt).targ1Pos(1);
            behView.targ1_y(t4) = dd(tt).targ1Pos(2);
            behView.targ2_x(t4) = dd(tt).targ2Pos(1);
            behView.targ2_y(t4) = dd(tt).targ2Pos(2);
            behView.t_targPres(t4) = dd(tt).targetsOn;
            behView.pulseTrial(t4) = dd(tt).pulseTrial;
            
            behView.t_pulseon(t4) = dd(tt).pulseOn;
            behView.t_pulseoff(t4) = dd(tt).pulseOff;
            behView.pulseSize(t4) = dd(tt).pulseSize;
            behView.pulseDur(t4) = dd(tt).pulseDuration;
            if isnan(dd(tt).pulseSize)
                behView.pulseCoh(t4) = dd(tt).sCoh;
            else
                behView.pulseCoh(t4) = dd(tt).pulseSize/1000 + dd(tt).sCoh;
            end
            behView.trialid(t4) = dd(tt).trialType;
            behView.rewDur(t4) = dd(tt).rewardSize;
        end
    end
end



end

