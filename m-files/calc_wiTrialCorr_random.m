
%% Project data onto first 46 "truly "random weight vectors (random unit vector) in each session (to generate 1034 unique pairs)

numPerms = 46;
proj_type = 'gauss';    % Method: Random weight vectors were generated from sampling from a Gaussian distribution
load(fullfile(saveLoc, ['mediation_randProj_' proj_type]), 'med_randProj', 't_med', 'tpts', 'randProj')

tic
clear perm_proj
for sess = 1 : 8
    
    disp(['Session ' num2str(sess)])
    % Load raw data for this session
    clear spk par dat_raw perm_proj_sess
    load(fullfile(saveLoc, ['spk_S' num2str(sess)]), 'spk', 'par')
    dat_raw = spk.SLnan;
    clear spk
    
    % behavioral indicators for all trials
    choice = par.beh.cho_trg';
    RT = par.beh.rt';
    coh = par.beh.sig_coh';
    
    % Project data onto the drections in state space defined by randProj{sess}
    for d = 1 : numPerms
        
        if rem(d, 25) == 0; disp(['Projecting Sess ' num2str(sess) ', d ' num2str(d)]); toc; end
        w_rand = randProj{sess}(d, :);
        
        clear dat_proj
        dat_proj = squeeze(mmx('mult', w_rand, dat_raw));
        
        perm_proj_sess(d, :, :) = dat_proj;
        
        if rem(d, 100) == 0
            disp(['Projecting Sess ' num2str(sess) ', d ' num2str(d)]); toc; 
            save(fullfile(saveLoc, ['perm_proj_sess' num2str(sess)]), 'perm_proj_sess', '-v7.3')
        end
    end
    disp('Saving...')
    save(fullfile(saveLoc, ['perm_proj_sess' num2str(sess)]), 'perm_proj_sess', '-v7.3')
    toc
end
disp('Finished projecting')


%% Calculate residuals of projected data with respect to mean activity per coherence
clear r_conv50b_rand
for sess = 1 : 8
    tic
    disp(['Generating residuals Sess ' num2str(sess) ', d ' num2str(d)]); toc; 
    
    clear perm_proj_sess r_conv50b_rand_sess 
    load(fullfile(saveLoc, ['perm_proj_sess' num2str(sess)]));
    disp('Done loading.'); toc
    
    par = par_sess{sess}; 
    
    sig_coh_sess = par.beh.sig_coh';    % Signed coherences 
    cohs = unique(sig_coh_sess);        % List of unique coherences
    
    numProj = size(perm_proj_sess, 1);
    for d = 1 : numProj
        
        if rem(d, 50) == 0; disp(['Creating residuals Sess ' num2str(sess) ', d ' num2str(d)]); toc; end
        % Create residuals
        dMat = squeeze(perm_proj_sess(d, :, :));
        if size(dMat, 1) > 0
            rMat = nan(size(dMat));
            for c = 1 : length(cohs)
                coh = cohs(c);
                trls = find(sig_coh_sess == coh);
                
                rMat(trls, :) = dMat(trls, :) - repmat(nanmean(dMat(trls, :), 1), length(trls), 1); % Compute residual within coherence
            end
            r_conv50b_rand_sess(d, :, :) = conv2(1, filtW50b, rMat, 'same'); % filtW50b: 50ms boxcar filter defined in general_settings
        end
    end
    disp('Saving...')
    save(fullfile(saveLoc, ['r_conv50b_rand_gauss_sess' num2str(sess)]), 'r_conv50b_rand_sess', '-v7.3')

end
clear perm_proj_sess r_conv50b_rand_sess


%% Compute within-trial correlation for all unique pairs between projectsions 
% for all sessions
tic
clear resCorr_rand
rtAdd = 0.1;
load(fullfile(saveLoc, 'par_sess'), 'par_sess')
time = par_sess{sess}.tsl';
for sess = 1 : 8
    
    disp(['Correlating Sess ' num2str(sess)]); toc;
    clear r_conv50b_rand_sess
    load(fullfile(saveLoc, ['r_conv50b_rand_gauss_sess' num2str(sess)]))
    disp('Finished loading'); toc
    
    rt = par_sess{sess}.beh.rt;
    dd = 0;
    for dim1 = 1 : size(r_conv50b_rand_sess, 1)
        
        for dim2 = 1 : size(r_conv50b_rand_sess, 1)
            if dim1 > dim2
                dd = dd + 1;
                if rem(dd, 10) == 0; disp(['Creating residuals Sess ' num2str(sess) ', d ' num2str(dd)]); toc; end
                for trl = 1 : size(r_conv50b_rand_sess, 2)
                    tpts = find(~isnan(squeeze(r_conv50b_rand_sess(dim1, trl, :))) & ~isnan(squeeze(r_conv50b_rand_sess(dim2, trl, :))) &...
                        time >= 0.2 & time <= rt(trl) - rtAdd);
                    
                    if length(tpts) > 5
                        tpts2 = tpts(1 : 50 : end);
                        if length(tpts2) >= 4
                            [R, P] = corrcoef(r_conv50b_rand_sess(dim1, trl, tpts2), r_conv50b_rand_sess(dim2, trl, tpts2), 'rows', 'complete');
                            resCorr_rand{sess}(dd, trl) = R(1, 2);
                        end
                    end
                end
            end
        end
    end
    save(fullfile(saveLoc, 'resCorr_rand_gauss'), 'resCorr_rand', '-v7.3')
end


%% Extract the most important values and sort by dimesnion order
clear allMeds_rand all_zrMean_rand rz_sess_rand
for dim1 = 1 : size(resCorr_rand{1}, 1)-1
    for sess =  1 : 8
        d1 = resCorr_rand{sess}(dim1, :);
        d1_n = d1(~isnan(d1)); r_mean = averageCorrelation(d1_n);  % This one
        all_zrMean_rand(dim1, sess) = r_mean; % This one
    end
    rz_sess_rand(dim1) = nanmean(all_zrMean_rand(dim1, :), 2); % This one
end

save(fullfile(saveLoc, 'resCorr_wiTrial_rand_gauss'), 'par_sess', 'filtW50b', 'rz_sess_rand', '-v7.3');
% save(fullfile(saveLoc, 'resCorr_wiTrial_rand'), 'par_sess', 'filtW50b', 'resCorr', 'r_conv50b_rand', 'rz_sess', '-v7.3');
notes = ['4 sessions'];
% load(fullfile(saveLoc, 'resCorr_wiTrial'), 'rz_sess', 'notes')

disp(['Random projections with rand gauss'])
disp(['Across trial correlation (mean+-std) = ' num2str(nanmean(rz_sess_rand)) ' +- ' num2str(nanstd(rz_sess_rand))])
% Random projections with rand gauss (46 projections and 1034 pairs)
% Across trial correlation (mean+-std) = 0.0003034 +- 0.047881


%% ==============================================================

%% Within-trial correlation for signals generated with permuted real dimensions
load(fullfile(saveLoc, 'sig_allSessions'), 'sigSL')
allW = {sigSL.w.ramp, sigSL.w.PC1, sigSL.w.TinC, sigSL.w.whatD, sigSL.w.whenC};
numPerms = 32;
clear sigSL 
cLen = 200;

tic
clear permute_proj_sess
for sess = 1 : 8
    
    % Load raw data for this session
    clear spk par dat_raw permute_proj_sess
    load(fullfile(saveLoc, ['spk_S' num2str(sess)]), 'spk', 'par')
    dat_raw = spk.SLnan;
    clear spk
    
    % Project data onto the drections in state space defined by randProj{sess}
    for dim = 1 : length(allW)
        
        w_orig = allW{dim}{sess};
        
        for d = 1 : numPerms
            
            if rem(d, 10) == 0; disp(['Projecting Sess ' num2str(sess) ', d ' num2str(d)]); toc; end
            p1 = randperm(length(w_orig));
            w_rand = w_orig(p1);
            
            permute_w{sess}(dim, d, :) = w_rand;
            
            clear dat_proj
            %dat_proj = squeeze(mmx('mult', w_rand, dat_raw));
            
            dat_proj = nan(size(dat_raw, 2), size(dat_raw, 3));
            for cc = 1 : ceil(size(dat_raw, 3)/cLen)
                tps = [((cc-1)*cLen + 1) : min(((cc-1)*cLen + cLen), size(dat_raw, 3))];
                dat_proj(:, tps) = squeeze(mmx('mult', w_rand, dat_raw(:, :, tps)));
            end
            
            permute_proj_sess{sess}(dim, d, :, :) = dat_proj;
        end
    end
    
    disp('Saving...')
    save(fullfile(saveLoc, ['permute_proj_sess' num2str(sess)]), 'permute_proj_sess', 'permute_w', '-v7.3')
    toc
end
disp('Finished projecting')


%% Calculate residuals of projected data with respect to mean activity per coherence
clear r_conv50b_rand_permute_sess
load(fullfile(saveLoc, 'par_sess'))
for sess = 1 : 8
    tic
    disp(['Generating residuals Sess ' num2str(sess) ', d ' num2str(d)]); toc;
    
    clear permute_proj_sess r_conv50b_rand_permute_sess
    load(fullfile(saveLoc, ['permute_proj_sess' num2str(sess)]));
    disp('Done loading.'); toc
    
    par = par_sess{sess};
    
    sig_coh_sess = par.beh.sig_coh';    % Signed coherences
    cohs = unique(sig_coh_sess);        % List of unique coherences
    
    numProj = size(permute_proj_sess{sess}, 2);
    numDim = size(permute_proj_sess{sess}, 1);
    for dim = 1 : numDim
        for d = 1 : numProj
            
            if rem(d, 10) == 0; disp(['Creating residuals Sess ' num2str(sess) ', dim ' num2str(dim) ', perm ' num2str(d)]); toc; end
            % Create residuals
            dMat = squeeze(permute_proj_sess{sess}(dim, d, :, :));
            if size(dMat, 1) > 0
                rMat = nan(size(dMat));
                for c = 1 : length(cohs)
                    coh = cohs(c);
                    trls = find(sig_coh_sess == coh);
                    
                    rMat(trls, :) = dMat(trls, :) - repmat(nanmean(dMat(trls, :), 1), length(trls), 1); % Compute residual within coherence
                end
                r_conv50b_rand_permute_sess(dim, d, :, :) = conv2(1, filtW50b, rMat, 'same'); % filtW50b: 50ms boxcar filter defined in general_settings
            end
        end
    end
    disp('Saving...')
    save(fullfile(saveLoc, ['r_conv50b_rand_permute_sess' num2str(sess)]), 'r_conv50b_rand_permute_sess', '-v7.3')
    
end
clear perm_proj_sess r_conv50b_rand_permute_sess


%% Compute within-trial correlation for all unique pairs between projectsions - shuffled vectors
% for all sessions
tic
clear resCorr_rand
rtAdd = 0.1;
load(fullfile(saveLoc, 'par_sess'), 'par_sess')
time = par_sess{sess}.tsl';
for sess = 1 : 8
    
    disp(['Correlating Sess ' num2str(sess)]); toc;
    clear r_conv50b_rand_sess
    load(fullfile(saveLoc, ['r_conv50b_rand_permute_sess' num2str(sess)]))
    disp('Finished loading'); toc
    
    rt = par_sess{sess}.beh.rt;
    
    numProj = size(r_conv50b_rand_permute_sess, 2);
    numDim = size(r_conv50b_rand_permute_sess, 1);
    numTrl = size(r_conv50b_rand_permute_sess, 3);
    
    
    for dim1 = 1 : numDim
        for dim2 = 1 : numDim
            if dim1 > dim2
                dd = 0;
                for d1 = 1 : numProj
                    for d2 = 1 : numProj
                        dd = dd + 1;
                        if rem(dd, 50) == 0; disp(['Correlating: Sess ' num2str(sess) ', d ' num2str(dd)]); toc; end
                        for trl = 1 : numTrl
                            tpts = find(~isnan(squeeze(r_conv50b_rand_permute_sess(dim1, d1, trl, :))) & ~isnan(squeeze(r_conv50b_rand_permute_sess(dim2, d2, trl, :))) &...
                                time >= 0.2 & time <= rt(trl) - rtAdd);
                            tpts2 = tpts(1 : 50 : end);
                            
                            if length(tpts2) >= 4
                                [R, P] = corrcoef(r_conv50b_rand_permute_sess(dim1, d1, trl, tpts2), r_conv50b_rand_permute_sess(dim2, d2, trl, tpts2), 'rows', 'complete');
                                resCorr_rand_permute{sess, dim1, dim2}(dd, trl) = R(1, 2);
                            end
                        end
                    end
                end
            end
        end
    end
    disp('Saving...'); toc
    save(fullfile(saveLoc, 'resCorr_rand_permute'), 'resCorr_rand_permute', '-v7.3')
    disp('Done saving'); toc
end











