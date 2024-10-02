function [spkZ, spkMapZ, spkViewZ, spkZnan, spkMapZnan, spkViewZnan] = normDat(spk, spkMap, spkView, par)

% Normalize firing rates based on mean and std of activity between 200 and
% 600 ms after motion onset during RDM discrimination trials
% Inputs
%   - par = basic information about the recording session and parameters
%           including 
%           the time vectors par.tsl (time aligned to motion onset in ms)
%           and par.trl (time aligned to motion onset in ms)
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
% Outputs
%  - [spkZ, spkMapZ, spkViewZ, spkZnan, spkMapZnan, spkViewZnan] =
%           normalized instances of the input variables (same format)


nTrials = size(spk.SL, 2);
nT = size(spk.SL, 3);
t_norm = [0.2 0.6]; 
t_incl = find(par.tsl >= t_norm(1) & par.tsl <= t_norm(2));
clear m s
for n = 1 : size(spk.SLnan, 1)
    
    mm = squeeze(spk.SLnan(n, :, t_incl));
    [~,mu,sigma] = zscore(mm(~isnan(mm)));
    
    s(n) = sigma;
    m(n) = mu;
end
m=m'; s=s';
% Exclusions of neurons that do not spike at all during reference period
goodN = find(s ~= 0);

% Normalize motion discrimination task data
mMat = repmat(m(goodN),[1 nTrials nT]) ;
sMat = repmat(s(goodN),[1 nTrials  nT]) ;
spkZnan = [];
spkZnan.m = m;
spkZnan.s = s;
spkZnan.goodN = goodN;
spkZnan.t_norm = t_norm;
spkZnan.SL = nan(size(spk.SLnan));
spkZnan.SL(goodN, :, :) =  (spk.SLnan(goodN, :, :) - mMat) ./ sMat ;
spkZnan.RL = nan(size(spk.RLnan));
spkZnan.RL(goodN, :, :) =  (spk.RLnan(goodN, :, :) - mMat) ./ sMat ;

spkZ = [];
spkZ.m = m;
spkZ.s = s;
spkZ.goodN = goodN;
spkZ.t_norm = t_norm;
spkZ.SL = nan(size(spk.SL));
spkZ.SL(goodN, :, :) =  (spk.SL(goodN, :, :) - mMat) ./ sMat ;
spkZ.RL = nan(size(spk.RL));
spkZ.RL(goodN, :, :) =  (spk.RL(goodN, :, :) - mMat) ./ sMat ;
clear spk

% Normalize mapping task data
nTrials = size(spkMap.SLnan, 2);
nT = size(spkMap.SLnan, 3);
mMat = repmat(m(goodN),[1 nTrials nT]) ;
sMat = repmat(s(goodN),[1 nTrials  nT]) ;
spkMapZnan = [];
spkMapZnan.m = m;
spkMapZnan.s = s;
spkMapZnan.goodN = goodN;
spkMapZnan.t_norm = t_norm;
spkMapZnan.SL = nan(size(spkMap.SLnan));
spkMapZnan.SL(goodN, :, :) =  (spkMap.SLnan(goodN, :, :) - mMat) ./ sMat ;
spkMapZnan.RL = nan(size(spkMap.RLnan));
spkMapZnan.RL(goodN, :, :) =  (spkMap.RLnan(goodN, :, :) - mMat) ./ sMat ;

spkMapZ = [];
spkMapZ.m = m;
spkMapZ.s = s;
spkMapZ.goodN = goodN;
spkMapZ.t_norm = t_norm;
spkMapZ.SL = nan(size(spkMap.SL));
spkMapZ.SL(goodN, :, :) =  (spkMap.SL(goodN, :, :) - mMat) ./ sMat ;
spkMapZ.RL = nan(size(spkMap.RL));
spkMapZ.RL(goodN, :, :) =  (spkMap.RL(goodN, :, :) - mMat) ./ sMat ;

% Normalize motion viewing task data
nTrials = size(spkView.SLnan, 2);
nT = size(spkView.SLnan, 3);
mMat = repmat(m(goodN),[1 nTrials nT]) ;
sMat = repmat(s(goodN),[1 nTrials  nT]) ;
spkViewZnan = [];
spkViewZnan.m = m;
spkViewZnan.s = s;
spkViewZnan.goodN = goodN;
spkViewZnan.t_norm = t_norm;
spkViewZnan.SL = nan(size(spkView.SLnan));
spkViewZnan.SL(goodN, :, :) =  (spkView.SLnan(goodN, :, :) - mMat) ./ sMat ;

spkViewZ = [];
spkViewZ.m = m;
spkViewZ.s = s;
spkViewZ.goodN = goodN;
spkViewZ.t_norm = t_norm;
spkViewZ.SL = nan(size(spkView.SL));
spkViewZ.SL(goodN, :, :) =  (spkView.SL(goodN, :, :) - mMat) ./ sMat ;

end