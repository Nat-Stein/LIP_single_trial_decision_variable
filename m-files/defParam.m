function [par] = defParam(sess, dat_path)
% loadParam: assign parameters for neuropixel sessions
%   Detailed explanation goes here

switch sess
    case 1
        par.date = '201211M';
        sess = 1;
        par.recordingName = '20201211_Mars_g0_t2_LIP';
        
    case 2
        par.date = '201030M';
        par.recordingName = '20201030_Mars_g0_t1_LIP';
%         par.rawPath = 'I:\Gabe\20201030';
    case 3
        par.date = '201110M';
        par.recordingName = '20201110_Mars_g0_t1_LIP';
%         par.rawPath = 'I:\Gabe\20201110';
    case 4
        par.date = '201116M';
        par.recordingName = '20201116_Mars_g0_t1_LIP';
%         par.rawPath = 'I:\Gabe\20201116';
    case 5
        par.date = '201208M';
        par.recordingName = '20201208_Mars_g0_t1_LIP';
%         par.rawPath = 'I:\Gabe\20201208';
    case 6
        par.date = '211011J';
        par.recordingName = '20211011_Jones_2_g0_t1_LIP';
%         par.rawPath = 'I:\Gabe\20211011';
    case 7
        par.date = '211015J';
        par.recordingName = '20211015_Jones_g0_t1_LIP'; 
%         par.rawPath = 'I:\Gabe\20211015'; 
    case 8
        par.date = '211020J';
        par.recordingName = '20211020_Jones_g0_t0_LIP'; 
%         par.rawPath = 'I:\Gabe\20211020'; 
    otherwise
        error('Session does not exist: %s', sess)
end
par.rawPath = fullfile(dat_path, par.recordingName(1:8));
par.datPath = dat_path;

par.matPath = fullfile(par.rawPath, 'mat-files'); % mat files will be saved here
par.saveFig = fullfile(par.rawPath, 'Figures'); % Single-unit figures will be saved here
par.figPath = fullfile(par.rawPath, 'Figures');  % Other figures files will be saved here
par.saveFigPop = fullfile(par.saveFig, 'Pop figures');  % Figures of population activity will be saved here
par.fileName = [par.recordingName '.mat'];

par.LIP_MIP = 1; % Get rid of this?

par.isihistbins = [0:0.001:2];
par.minISI = 0.001;

par.eplims = {[-400 4000] / 1000, [-4000 400] / 1000, [-300 500] / 1000};
par.dt = 0.001;
par.tsl = par.eplims{1}(1) : par.dt : par.eplims{1}(2);
par.trl = par.eplims{2}(1) : par.dt : par.eplims{2}(2);
par.trg = par.eplims{3}(1) : par.dt : par.eplims{3}(2);
par.tfp = par.eplims{3}(1) : par.dt : par.eplims{3}(2);

par.cutSL = 0.1; % time cut before saccade in stimulus-locked epochs [s]
par.cutRL = 0.14; % time cut after stimulus onset in response locked epochs [s]
par.cutFPoff = 0;
par.cutRLMap = 0.3;

end


