function [dd, par] = loadSess(par)
%loadSess load data of single session for preprocessing


matName = fullfile(par.rawPath,[par.recordingName]);
load(matName);
dd = d; clear d

load(fullfile(par.rawPath, 'unitIdxLIP_Tin'))
par.unit_Tin = unitIdxLIP_Tin;
load(fullfile(par.rawPath, 'unitIdxLIP_Tout'))
par.unit_Tout = unitIdxLIP_Tout;

% RENAME AND FILL IN NUMBERS
% par.unit_MinC = []; % unitIdx_DinRFcC
% par.unit_MinI = []; % unitIdx_DinRFiC
% par.unit_MinN = par.unitIdx_DinRFnC;
end

