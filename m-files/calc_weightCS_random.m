% Calculate cosine similarity between random projections in state space

load(fullfile(saveLoc, 'mediation_randProj_gauss'), 'randProj')

numPerm = size(randProj{1},1);
cosine_sim_rand = nan(8,numPerm);
for sess = 1 : 8
    for n = 1 : numPerm - 1
        w1 = randProj{sess}(n, :);
        w2 = randProj{sess}(n+1, :);
        
        xy   = dot(w1,w2);
        nx   = norm(w1);
        ny   = norm(w2);
        nxny = nx*ny;
        Cs   = xy/nxny;
        cosine_sim_rand(sess, n) = Cs;
    end
    disp(['CS Sess ' num2str(sess) ': mean+-std= ' num2str(nanmean(cosine_sim_rand(sess, :))) '+-' num2str(nanstd(cosine_sim_rand(sess, :)))])
    
end

save(fullfile(saveLoc, 'cosine_sim_rand'), 'cosine_sim_rand')

disp(['CS across sessions: mean+-std(across permutations)= ' num2str(nanmean(nanmean(cosine_sim_rand, 2), 1)) '+-' num2str(nanstd(nanmean(cosine_sim_rand, 1)))])
% Uniform distribution: 
% CS across sessions: mean+-sem(across permutations)= 0.0054002+-0.0010701
% Gaussian distribution:
% CS Sess 1: mean+-std= 0.0003018+-0.072152
% CS Sess 2: mean+-std= -0.0043596+-0.10211
% CS Sess 3: mean+-std= -0.0034649+-0.13839
% CS Sess 4: mean+-std= 0.00030166+-0.085141
% CS Sess 5: mean+-std= 0.0013611+-0.095582
% CS Sess 6: mean+-std= -0.0017052+-0.081527
% CS Sess 7: mean+-std= 0.00061816+-0.078316
% CS Sess 8: mean+-std= 0.0029+-0.069427
% CS across sessions: mean+-std(across permutations)= -0.00050588+-0.032598


% Create null distributions by shuffling weights of existing vectors
allW = {sigSL.w.ramp, sigSL.w.PC1, sigSL.w.TinC, sigSL.w.whatD, sigSL.w.whenC};

cosine_sim_randperm = nan(8, 5, 5, numPerm);
for sess = 1 : 8
    for d1 = 1 : length(allW)
        w1_orig = allW{d1}{1, sess};
        numN = length(w1_orig);
        for d2 = 1 : length(allW)
            w2_orig = allW{d2}{1, sess};
            if d1 ~= d2
                for n = 1 : numPerm - 1
                    p1 = randperm(numN);
                    w1 = w1_orig(p1);
                    
                    p2 = randperm(numN);
                    w2 = w2_orig(p2);
                    
                    xy   = dot(w1,w2);
                    nx   = norm(w1);
                    ny   = norm(w2);
                    nxny = nx*ny;
                    Cs   = xy/nxny;
                    cosine_sim_randperm(sess, d1, d2, n) = Cs;
                end
            end
        end
    end
    disp(['CS Sess ' num2str(sess) ': mean+-std= ' num2str(nanmean(cosine_sim_rand(sess, :))) '+-' num2str(nanstd(cosine_sim_rand(sess, :)))])
    
end
save(fullfile(saveLoc, 'cosine_sim_randperm'), 'cosine_sim_randperm') % cosine_sim_randperm(sess, d1, d2, :)


for sess = 1 : 8
    for d1 = 1 : length(allW)
        for d2 = 1 : length(allW)
            
            disp(['CS Sess ' num2str(sess) ': mean+-std= ' num2str(nanmean(cosine_sim_randperm(sess, d1, d2, :))) '+-'...
                num2str(nanstd(cosine_sim_randperm(sess, d1, d2, :)))])
            
        end
    end
end

figure; imagesc(squeeze(nanmean(nanmean(cosine_sim_randperm(:, :, :, :), 1), 4)))




