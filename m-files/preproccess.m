% Preprocess neuropixel data
% Extract and reformat behavioral data and spikes for each session
% Output structures:
% beh, behMap, behView: behavioral indicators by trial in RDM, mapping
% and motion viewing tasks
% spk, spkMap, spkView: neural data by trial in RDM, mapping
% and motion viewing tasks


%% Load data of individual session
for sess = 1 : 8
    
    disp(['Pre-processing session ' num2str(sess)])
    clear par dd beh behMap behView spk spkMap spkView
    clear spkZ spkMapZ spkViewZ spkZnan spkMapZnan spkViewZnan
    par = defParam(sess, dat_path);
    
    %% Load data of single session
    [dd, par] = loadSess(par);
    
    %% Re-format behavioral data
    [beh, behMap, behView] = reformBehav(par, dd);
    
    par.beh = beh;
    par.behMap = behMap;
    par.behView = behView;
    
    disp(['Saving session ' num2str(sess) ' - behaviour'])
    save(fullfile(saveLoc, ['beh_S' num2str(sess)]), 'beh', 'behMap', 'behView', 'par')
    %% Reformat spiking data
    [spk, spkMap, spkView, par] = reformDat(par, dd);
    
    disp(['Saving session ' num2str(sess)  ' - neural activity'])
    save(fullfile(saveLoc, ['spk_S' num2str(sess)]), 'spk', 'spkMap', 'spkView', 'par', '-v7.3')
    
    %% Normalize spiking data
    disp(['Computing normalization: session ' num2str(sess)])
    [spkZ, spkMapZ, spkViewZ, spkZnan, spkMapZnan, spkViewZnan] = normDat(spk, spkMap, spkView, par);
    disp(['Saving session ' num2str(sess) ' - normalized neural activity'])
    save(fullfile(saveLoc, ['spkZ_S' num2str(sess)]),...
        'spkZ', 'spkMapZ', 'spkViewZ','spkZnan', 'spkMapZnan', 'spkViewZnan', 'par', '-v7.3')
    
end

