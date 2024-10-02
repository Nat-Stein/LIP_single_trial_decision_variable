% comb_par

load(fullfile(saveLoc, ['spk_S' num2str(1)]), 'par')
% Create shorter time vector
tsl_s = [-0.2 1]; tpts = find(par.tsl >= tsl_s(1) & par.tsl <= tsl_s(2));
trl_s = -1 * flip(tsl_s); tpts_r = find(par.trl >= trl_s(1) & par.trl <= trl_s(2));
tsl = par.tsl(tpts);
trl = par.trl(tpts_r);

global par_sess
for sess = 1 : 8
    % Load single session parameter structure
    clear par
    load(fullfile(saveLoc, ['spk_S' num2str(sess)]), 'par')
    
    par.datPath = dat_path;
    par.saveLoc = saveLoc;
    
    % Load predefined Min neurons - code still needs to be adapted
    %loadPredefinedMins
    
    nCell = par.numCells; % Is this included now?
    %% Derive weight vectors for
    for d = [1 2 7 8]
        dim = dNames{d};
        w = zeros(1, nCell);
        switch d
            case 1 % TinC
                w(par.unit_Tin) = 1/length(par.unit_Tin); 
            case 2 % TinI
                w(par.unit_Tout) = 1/length(par.unit_Tout); 
            case 7 % MinC
                par.unit_MinC = find(w > 0); % Placeholder
            case 8 % MinI
                %w = w_minI; % Placeholder
                par.unit_MinI = find(w > 0); % Placeholder
        end
        ff = fullfile(saveLoc,...
            [dNames{d} '_w_S' num2str(sess)]);
        save(ff, 'w')
        if size(w, 1) > size(w, 2); w = w'; save(ff, 'w'); end
    end
    par_sess{sess} = par;
end
save(fullfile(saveLoc, 'par_sess'), 'par_sess')