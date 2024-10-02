function [S,P]  = simRace2DDM(t, varargin)
% returns stucture S with fields
% .decTime (decision timer)
% .choice (0 or 1)
% .dvTin (dv of the positive choice)

% test 
% t = (0:0.001:1); [X,P] = simRace2DDM(t);
% simRace2DDM(t,'seed',P.seed);

% varargin
pSet = inputParser;
addParameter(pSet,'coh',0); % -1..1 scale
addParameter(pSet,'kappa',8);
addParameter(pSet,'B',0.9);
addParameter(pSet,'BloFactor',1.5); % reflecting lower bound in units of -B
addParameter(pSet,'UrgLin',0); % slope of urgency (linear)
addParameter(pSet,'r',-.7, @(x) x<=0 && x>=-1); % correlation
addParameter(pSet,'PlotTraces',false,@islogical);
addParameter(pSet,'seed',nan);
addParameter(pSet,'untermTval',nan)
% addParameter(pSet,'Q', nan, @ismatrix);
% option to freeze Rand num seed

parse(pSet,varargin{:});
P = pSet.Results;

% check
if t(1) ~= 0 || t(2)~=0.001
    error('t must start at zero and dt must be 0.001');
else
    dt = 0.001;
end
B = P.B;
Bref = -B * P.BloFactor; % reflecting lower bound
% if isnan(P.Blo)
%     Bref = -P.B; % reflecting lower bound
% else
%     Bref = P.Blo;
% end

% variables
Nt = length(t);
Q = sqrtm([1 P.r; P.r 1]); 
    
% dvPos = nan(Nt,1);
dvPos = nan(size(t));
dvNeg = nan(size(t));

if isstruct(P.seed)
    rng(P.seed);
end

% deterministic components
d_drift = P.kappa*P.coh*dt; % drift per dt
d_urg = P.UrgLin*dt; % urgeny per dt

% generate the anticorrelated momentary evidence for this trial
w = randn(Nt,2) * Q * sqrt(dt); % Columns are pos and neg race 

% integrate both races a step at a time
i=1; isTerm = false;
dvPos(1) = 0; dvNeg(1)=0;
choice = nan;
while i<Nt && ~isTerm
    i = i+1; %  start at i=2
    dvPos(i) = dvPos(i-1) +  w(i,1) + d_drift + d_urg;
    if dvPos(i) <= Bref
        dvPos(i) = Bref;
    elseif dvPos(i) >= B
        isTerm = true;
        S.decTime = t(i);
        S.choice = 1;
    end
    dvNeg(i) = dvNeg(i-1) + w(i,2) - d_drift + d_urg;
    if dvNeg(i) <= Bref
        dvNeg(i) = Bref;
    elseif dvNeg(i) >= B
        isTerm = true;
        S.decTime = t(i);
        S.choice = 0;
    end
end
S.term = isTerm;
if ~isTerm    
    if dvPos(i) > dvNeg(i)
        S.choice = 1;
    else
        S.choice=0;
    end
    S.decTime = P.untermTval; % nan is default. Consider max of t + 0.2?
end
S.dvTin = dvPos;
S.dvTout = dvNeg;
if P.PlotTraces
    plot(t,S.dvTin,"k-"), hold on
    plot(t,S.dvTout,"k--"), hold off
end
P.seed = rng; % useful for repeating identical random numbers



