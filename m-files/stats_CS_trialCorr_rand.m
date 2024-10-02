
load(fullfile(saveLoc, 'cosine_sim'), 'cosine_sim')
load(fullfile(saveLoc, 'resCorr_wiTrial'), 'all_zrMean', 'resCorr');
load(fullfile(saveLoc, 'cosine_sim_randperm'), 'cosine_sim_randperm') % cosine_sim_randperm(sess, d1, d2, :)
load(fullfile(saveLoc, 'resCorr_wiTrial_rand_gauss'), 'par_sess', 'filtW50b', 'rz_sess_rand')
load(fullfile(saveLoc, 'cosine_sim_rand'), 'cosine_sim_rand')
load(fullfile(saveLoc, 'resCorr_rand_permute'), 'resCorr_rand_permute')




%% Cosine similarities stats: data vs random (permuted)

load(fullfile(saveLoc, 'cosine_sim'), 'cosine_sim')

allpVals = [];
for d1 = 1 : length(dimLabelsCS)
    for d2 = 1 : length(dimLabelsCS)
        if d1 < d2
            for sess = 1 : 8
                cs_data = squeeze(cosine_sim(sess, d1, d2));
                cs_randperm = squeeze(cosine_sim_randperm(sess, d1, d2, :));
                [h, p, ci, tStat] = ttest(cs_randperm, cs_data);
                disp(['S' num2str(sess) ': ' dimLabelsCS{d1} ' vs ' dimLabelsCS{d2} ': p=' num2str(p)])
                %disp(['mean +- sem = ' num2str(nanmean(cs_randperm)) ' +- ' num2str(nanstd(cs_randperm))])
                allpVals(d1, d2, sess) = p;
            end
        end
    end
end
disp(['All p < ' num2str(max(max(max(allpVals))))]); % All p < 1.1793e-157

% Output
% S1: ramp vs PC1: p=0
% S2: ramp vs PC1: p=0
% S3: ramp vs PC1: p=0
% S4: ramp vs PC1: p=0
% S5: ramp vs PC1: p=0
% S6: ramp vs PC1: p=0
% S7: ramp vs PC1: p=0
% S8: ramp vs PC1: p=0
% S1: ramp vs T^{con}_{in}: p=0
% S2: ramp vs T^{con}_{in}: p=5.6124e-240
% S3: ramp vs T^{con}_{in}: p=0
% S4: ramp vs T^{con}_{in}: p=0
% S5: ramp vs T^{con}_{in}: p=0
% S6: ramp vs T^{con}_{in}: p=0
% S7: ramp vs T^{con}_{in}: p=0
% S8: ramp vs T^{con}_{in}: p=0
% S1: ramp vs What-decoder: p=0
% S2: ramp vs What-decoder: p=0
% S3: ramp vs What-decoder: p=0
% S4: ramp vs What-decoder: p=0
% S5: ramp vs What-decoder: p=0
% S6: ramp vs What-decoder: p=0
% S7: ramp vs What-decoder: p=0
% S8: ramp vs What-decoder: p=0
% S1: ramp vs When-decoder: p=0
% S2: ramp vs When-decoder: p=0
% S3: ramp vs When-decoder: p=0
% S4: ramp vs When-decoder: p=0
% S5: ramp vs When-decoder: p=0
% S6: ramp vs When-decoder: p=0
% S7: ramp vs When-decoder: p=0
% S8: ramp vs When-decoder: p=0
% S1: PC1 vs T^{con}_{in}: p=0
% S2: PC1 vs T^{con}_{in}: p=3.3723e-220
% S3: PC1 vs T^{con}_{in}: p=0
% S4: PC1 vs T^{con}_{in}: p=0
% S5: PC1 vs T^{con}_{in}: p=0
% S6: PC1 vs T^{con}_{in}: p=0
% S7: PC1 vs T^{con}_{in}: p=0
% S8: PC1 vs T^{con}_{in}: p=0
% S1: PC1 vs What-decoder: p=0
% S2: PC1 vs What-decoder: p=0
% S3: PC1 vs What-decoder: p=0
% S4: PC1 vs What-decoder: p=3.7336e-183
% S5: PC1 vs What-decoder: p=4.578e-320
% S6: PC1 vs What-decoder: p=0
% S7: PC1 vs What-decoder: p=0
% S8: PC1 vs What-decoder: p=0
% S1: PC1 vs When-decoder: p=0
% S2: PC1 vs When-decoder: p=0
% S3: PC1 vs When-decoder: p=0
% S4: PC1 vs When-decoder: p=0
% S5: PC1 vs When-decoder: p=0
% S6: PC1 vs When-decoder: p=0
% S7: PC1 vs When-decoder: p=0
% S8: PC1 vs When-decoder: p=0
% S1: T^{con}_{in} vs What-decoder: p=0
% S2: T^{con}_{in} vs What-decoder: p=3.312e-312
% S3: T^{con}_{in} vs What-decoder: p=1.7357e-269
% S4: T^{con}_{in} vs What-decoder: p=0
% S5: T^{con}_{in} vs What-decoder: p=0
% S6: T^{con}_{in} vs What-decoder: p=0
% S7: T^{con}_{in} vs What-decoder: p=0
% S8: T^{con}_{in} vs What-decoder: p=6.7573e-290
% S1: T^{con}_{in} vs When-decoder: p=0
% S2: T^{con}_{in} vs When-decoder: p=1.5307e-228
% S3: T^{con}_{in} vs When-decoder: p=1.1793e-157
% S4: T^{con}_{in} vs When-decoder: p=0
% S5: T^{con}_{in} vs When-decoder: p=0
% S6: T^{con}_{in} vs When-decoder: p=0
% S7: T^{con}_{in} vs When-decoder: p=0
% S8: T^{con}_{in} vs When-decoder: p=0
% S1: What-decoder vs When-decoder: p=0
% S2: What-decoder vs When-decoder: p=0
% S3: What-decoder vs When-decoder: p=0
% S4: What-decoder vs When-decoder: p=0
% S5: What-decoder vs When-decoder: p=0
% S6: What-decoder vs When-decoder: p=0
% S7: What-decoder vs When-decoder: p=0
% S8: What-decoder vs When-decoder: p=0
% All p < 1.1793e-157


%% Cosine similarities stats: data vs random (random unit vectors)

allpVals = [];
for d1 = 1 : length(dimLabelsCS)
    for d2 = 1 : length(dimLabelsCS)
        if d1 < d2
            for sess = 1 : 8
                cs_data = squeeze(cosine_sim(sess, d1, d2));
                cs_rand = squeeze(cosine_sim_rand(sess, :));
                [h, p, ci, tStat] = ttest(cs_rand, cs_data);
                disp(['S' num2str(sess) ': ' dimLabelsCS{d1} ' vs ' dimLabelsCS{d2} ': p=' num2str(p)])
                %disp(['mean +- sem = ' num2str(nanmean(cs_randperm)) ' +- ' num2str(nanstd(cs_randperm))])
                allpVals(d1, d2, sess) = p;
            end
        end
    end
end
disp(['All p < ' num2str(max(max(max(allpVals))))]); % All p < 2.0166e-64

% Output
% S1: ramp vs PC1: p=0
% S2: ramp vs PC1: p=0
% S3: ramp vs PC1: p=0
% S4: ramp vs PC1: p=0
% S5: ramp vs PC1: p=0
% S6: ramp vs PC1: p=0
% S7: ramp vs PC1: p=0
% S8: ramp vs PC1: p=0
% S1: ramp vs T^{con}_{in}: p=0
% S2: ramp vs T^{con}_{in}: p=8.6739e-183
% S3: ramp vs T^{con}_{in}: p=0
% S4: ramp vs T^{con}_{in}: p=0
% S5: ramp vs T^{con}_{in}: p=0
% S6: ramp vs T^{con}_{in}: p=0
% S7: ramp vs T^{con}_{in}: p=0
% S8: ramp vs T^{con}_{in}: p=0
% S1: ramp vs What-decoder: p=0
% S2: ramp vs What-decoder: p=0
% S3: ramp vs What-decoder: p=0
% S4: ramp vs What-decoder: p=0
% S5: ramp vs What-decoder: p=0
% S6: ramp vs What-decoder: p=0
% S7: ramp vs What-decoder: p=0
% S8: ramp vs What-decoder: p=0
% S1: ramp vs When-decoder: p=0
% S2: ramp vs When-decoder: p=0
% S3: ramp vs When-decoder: p=0
% S4: ramp vs When-decoder: p=0
% S5: ramp vs When-decoder: p=0
% S6: ramp vs When-decoder: p=0
% S7: ramp vs When-decoder: p=0
% S8: ramp vs When-decoder: p=0
% S1: PC1 vs T^{con}_{in}: p=0
% S2: PC1 vs T^{con}_{in}: p=2.0166e-64
% S3: PC1 vs T^{con}_{in}: p=0
% S4: PC1 vs T^{con}_{in}: p=0
% S5: PC1 vs T^{con}_{in}: p=0
% S6: PC1 vs T^{con}_{in}: p=0
% S7: PC1 vs T^{con}_{in}: p=0
% S8: PC1 vs T^{con}_{in}: p=0
% S1: PC1 vs What-decoder: p=0
% S2: PC1 vs What-decoder: p=0
% S3: PC1 vs What-decoder: p=0
% S4: PC1 vs What-decoder: p=1.8075e-132
% S5: PC1 vs What-decoder: p=1.6008e-321
% S6: PC1 vs What-decoder: p=0
% S7: PC1 vs What-decoder: p=0
% S8: PC1 vs What-decoder: p=0
% S1: PC1 vs When-decoder: p=0
% S2: PC1 vs When-decoder: p=0
% S3: PC1 vs When-decoder: p=0
% S4: PC1 vs When-decoder: p=0
% S5: PC1 vs When-decoder: p=0
% S6: PC1 vs When-decoder: p=0
% S7: PC1 vs When-decoder: p=0
% S8: PC1 vs When-decoder: p=0
% S1: T^{con}_{in} vs What-decoder: p=0
% S2: T^{con}_{in} vs What-decoder: p=8.8235e-306
% S3: T^{con}_{in} vs What-decoder: p=0
% S4: T^{con}_{in} vs What-decoder: p=0
% S5: T^{con}_{in} vs What-decoder: p=0
% S6: T^{con}_{in} vs What-decoder: p=0
% S7: T^{con}_{in} vs What-decoder: p=0
% S8: T^{con}_{in} vs What-decoder: p=2.2702e-293
% S1: T^{con}_{in} vs When-decoder: p=0
% S2: T^{con}_{in} vs When-decoder: p=1.1193e-290
% S3: T^{con}_{in} vs When-decoder: p=8.8269e-252
% S4: T^{con}_{in} vs When-decoder: p=0
% S5: T^{con}_{in} vs When-decoder: p=0
% S6: T^{con}_{in} vs When-decoder: p=0
% S7: T^{con}_{in} vs When-decoder: p=0
% S8: T^{con}_{in} vs When-decoder: p=0
% S1: What-decoder vs When-decoder: p=0
% S2: What-decoder vs When-decoder: p=0
% S3: What-decoder vs When-decoder: p=0
% S4: What-decoder vs When-decoder: p=0
% S5: What-decoder vs When-decoder: p=0
% S6: What-decoder vs When-decoder: p=0
% S7: What-decoder vs When-decoder: p=0
% S8: What-decoder vs When-decoder: p=0
% All p < 2.0166e-64

%% Within-trial correlation stats: data vs permuted weights


allPs = [];
for dim1 = 1 : size(resCorr, 1)
    for dim2 = 1 : size(resCorr, 2)
        if dim1 == dim2
            rz_sess(dim1, dim2) = 1;
        elseif dim1 > dim2
            for sess =  1 : 8
                rvals = resCorr{dim1, dim2, sess};
                
                rvals_data = rvals(~isnan(rvals));
                [r_mean_data, sem_corr_data, z_scores_data, p_data] = averageCorrelation(rvals_data);
                
                rvals_rand = resCorr_rand_permute{sess, dim1, dim2};
                rvals_rand = rvals_rand(:);
                rvals_rand = rvals_rand(~isnan(rvals_rand));
                [r_mean_rand, sem_corr_rand, z_scores_rand, p_rand] = averageCorrelation(rvals_rand');
                
                [h,p,ci,stats] = ttest2(z_scores_data, z_scores_rand);
                
                disp([dimLabelsCS{dim1} ' vs ' dimLabelsCS{dim2} ': p=' num2str(p)])
                allPs(sess, dim1, dim2) = p;
            end
            %rz_sess(dim1, dim2) = nanmean(all_zrMean{dim1, dim2}); % This one
        end
    end
end
disp(['All p <= ' num2str(max(max(max(allPs))))])

% Output
% PC1 vs ramp: p=0
% PC1 vs ramp: p=0
% PC1 vs ramp: p=0
% PC1 vs ramp: p=0
% PC1 vs ramp: p=0
% PC1 vs ramp: p=0
% PC1 vs ramp: p=0
% PC1 vs ramp: p=0
% Tin^{con}_{in} vs ramp: p=0
% Tin^{con}_{in} vs ramp: p=1.007e-39
% Tin^{con}_{in} vs ramp: p=0
% Tin^{con}_{in} vs ramp: p=7.012e-133
% Tin^{con}_{in} vs ramp: p=0
% Tin^{con}_{in} vs ramp: p=0
% Tin^{con}_{in} vs ramp: p=8.1878e-179
% Tin^{con}_{in} vs ramp: p=7.5916e-266
% Tin^{con}_{in} vs PC1: p=0
% Tin^{con}_{in} vs PC1: p=4.0399e-74
% Tin^{con}_{in} vs PC1: p=0
% Tin^{con}_{in} vs PC1: p=3.3702e-95
% Tin^{con}_{in} vs PC1: p=0
% Tin^{con}_{in} vs PC1: p=0
% Tin^{con}_{in} vs PC1: p=0
% Tin^{con}_{in} vs PC1: p=0
% What-decoder vs ramp: p=0
% What-decoder vs ramp: p=0
% What-decoder vs ramp: p=0
% What-decoder vs ramp: p=5.8306e-238
% What-decoder vs ramp: p=1.0883e-244
% What-decoder vs ramp: p=0
% What-decoder vs ramp: p=0
% What-decoder vs ramp: p=0
% What-decoder vs PC1: p=0
% What-decoder vs PC1: p=7.1403e-162
% What-decoder vs PC1: p=0
% What-decoder vs PC1: p=3.948e-100
% What-decoder vs PC1: p=1.3726e-176
% What-decoder vs PC1: p=0
% What-decoder vs PC1: p=0
% What-decoder vs PC1: p=0
% What-decoder vs Tin^{con}_{in}: p=0
% What-decoder vs Tin^{con}_{in}: p=9.1743e-130
% What-decoder vs Tin^{con}_{in}: p=0
% What-decoder vs Tin^{con}_{in}: p=1.5295e-65
% What-decoder vs Tin^{con}_{in}: p=8.9922e-123
% What-decoder vs Tin^{con}_{in}: p=0
% What-decoder vs Tin^{con}_{in}: p=0
% What-decoder vs Tin^{con}_{in}: p=0
% When-decoder vs ramp: p=0
% When-decoder vs ramp: p=0
% When-decoder vs ramp: p=0
% When-decoder vs ramp: p=0
% When-decoder vs ramp: p=0
% When-decoder vs ramp: p=0
% When-decoder vs ramp: p=0
% When-decoder vs ramp: p=0
% When-decoder vs PC1: p=0
% When-decoder vs PC1: p=0
% When-decoder vs PC1: p=0
% When-decoder vs PC1: p=2.9425e-96
% When-decoder vs PC1: p=0
% When-decoder vs PC1: p=0
% When-decoder vs PC1: p=0
% When-decoder vs PC1: p=0
% When-decoder vs Tin^{con}_{in}: p=0
% When-decoder vs Tin^{con}_{in}: p=1.302e-87
% When-decoder vs Tin^{con}_{in}: p=0
% When-decoder vs Tin^{con}_{in}: p=3.3006e-40
% When-decoder vs Tin^{con}_{in}: p=0
% When-decoder vs Tin^{con}_{in}: p=0
% When-decoder vs Tin^{con}_{in}: p=0
% When-decoder vs Tin^{con}_{in}: p=0
% When-decoder vs What-decoder: p=0
% When-decoder vs What-decoder: p=0
% When-decoder vs What-decoder: p=0
% When-decoder vs What-decoder: p=5.6062e-204
% When-decoder vs What-decoder: p=0
% When-decoder vs What-decoder: p=0
% When-decoder vs What-decoder: p=0
% When-decoder vs What-decoder: p=0
% All p <= 1.007e-39


%% Within-trial correlation stats: data vs random projections

allPs = [];
for dim1 = 1 : size(resCorr, 1)
    for dim2 = 1 : size(resCorr, 2)
        if dim1 == dim2
            rz_sess(dim1, dim2) = 1;
        elseif dim1 > dim2
            
            rvals_rand = rz_sess_rand;
            rvals_rand = rvals_rand(~isnan(rvals_rand));
            [r_mean_rand, sem_corr_rand, z_scores_rand, p_rand] = averageCorrelation(rvals_rand);
            
            for sess =  1 : 8
                rvals = resCorr{dim1, dim2, sess};
                
                rvals_data = rvals(~isnan(rvals));
                [r_mean_data, sem_corr_data, z_scores_data, p_data] = averageCorrelation(rvals_data);
                
                
                
                [h,p,ci,stats] = ttest2(z_scores_data, z_scores_rand);
                
                % Z-test
                % difference in means
                % sigma sqrt(sum(variances))
                
                % std_corr_rand = sem_corr_rand * sqrt(length(rvals_rand));

                %[h, p] = ztest(0, avg_corr_data, sem_corr_data);
                
                % Same as:
                % dataOut = r2z(resCorr{dim1, dim2, sess, 1});
                % z_mean = nanmean(dataOut);
                % r_mean = tanh(z_mean);
                %all_zrMean{dim1, dim2}(sess) = r_mean; % This one
                disp([dimLabelsCS{dim1} ' vs ' dimLabelsCS{dim2} ': p=' num2str(p)])
                allPs(sess, dim1, dim2) = p;
            end
            %rz_sess(dim1, dim2) = nanmean(all_zrMean{dim1, dim2}); % This one
        end
    end
end
disp(['All p <= ' num2str(max(max(max(allPs))))])

% Output
% PC1 vs ramp: p=0
% PC1 vs ramp: p=0
% PC1 vs ramp: p=0
% PC1 vs ramp: p=0
% PC1 vs ramp: p=0
% PC1 vs ramp: p=0
% PC1 vs ramp: p=2.6738e-200
% PC1 vs ramp: p=0
% Tin^{con}_{in} vs ramp: p=0
% Tin^{con}_{in} vs ramp: p=1.6083e-19
% Tin^{con}_{in} vs ramp: p=0
% Tin^{con}_{in} vs ramp: p=2.0922e-70
% Tin^{con}_{in} vs ramp: p=0
% Tin^{con}_{in} vs ramp: p=3.3483e-308
% Tin^{con}_{in} vs ramp: p=9.5144e-113
% Tin^{con}_{in} vs ramp: p=1.2787e-91
% Tin^{con}_{in} vs PC1: p=0
% Tin^{con}_{in} vs PC1: p=3.837e-25
% Tin^{con}_{in} vs PC1: p=1.1027e-298
% Tin^{con}_{in} vs PC1: p=1.0913e-33
% Tin^{con}_{in} vs PC1: p=0
% Tin^{con}_{in} vs PC1: p=0
% Tin^{con}_{in} vs PC1: p=0
% Tin^{con}_{in} vs PC1: p=0
% What-decoder vs ramp: p=5.7545e-297
% What-decoder vs ramp: p=7.4578e-281
% What-decoder vs ramp: p=0
% What-decoder vs ramp: p=4.6513e-122
% What-decoder vs ramp: p=4.62e-113
% What-decoder vs ramp: p=0
% What-decoder vs ramp: p=1.6032e-150
% What-decoder vs ramp: p=0
% What-decoder vs PC1: p=0
% What-decoder vs PC1: p=2.3811e-91
% What-decoder vs PC1: p=1.3154e-275
% What-decoder vs PC1: p=1.1139e-44
% What-decoder vs PC1: p=2.137e-83
% What-decoder vs PC1: p=0
% What-decoder vs PC1: p=0
% What-decoder vs PC1: p=0
% What-decoder vs Tin^{con}_{in}: p=0
% What-decoder vs Tin^{con}_{in}: p=1.5326e-76
% What-decoder vs Tin^{con}_{in}: p=1.5877e-245
% What-decoder vs Tin^{con}_{in}: p=1.5255e-53
% What-decoder vs Tin^{con}_{in}: p=1.0813e-80
% What-decoder vs Tin^{con}_{in}: p=0
% What-decoder vs Tin^{con}_{in}: p=0
% What-decoder vs Tin^{con}_{in}: p=2.8924e-147
% When-decoder vs ramp: p=0
% When-decoder vs ramp: p=0
% When-decoder vs ramp: p=1.1858e-321
% When-decoder vs ramp: p=1.5743e-225
% When-decoder vs ramp: p=0
% When-decoder vs ramp: p=0
% When-decoder vs ramp: p=0
% When-decoder vs ramp: p=0
% When-decoder vs PC1: p=0
% When-decoder vs PC1: p=7.4585e-172
% When-decoder vs PC1: p=0
% When-decoder vs PC1: p=2.1412e-49
% When-decoder vs PC1: p=0
% When-decoder vs PC1: p=0
% When-decoder vs PC1: p=0
% When-decoder vs PC1: p=0
% When-decoder vs Tin^{con}_{in}: p=0
% When-decoder vs Tin^{con}_{in}: p=7.7488e-63
% When-decoder vs Tin^{con}_{in}: p=3.5064e-316
% When-decoder vs Tin^{con}_{in}: p=3.2259e-36
% When-decoder vs Tin^{con}_{in}: p=0
% When-decoder vs Tin^{con}_{in}: p=0
% When-decoder vs Tin^{con}_{in}: p=7.8556e-314
% When-decoder vs Tin^{con}_{in}: p=2.8624e-152
% When-decoder vs What-decoder: p=0
% When-decoder vs What-decoder: p=0
% When-decoder vs What-decoder: p=1.0066e-307
% When-decoder vs What-decoder: p=5.0483e-123
% When-decoder vs What-decoder: p=8.6297e-207
% When-decoder vs What-decoder: p=0
% When-decoder vs What-decoder: p=0
% When-decoder vs What-decoder: p=0
% All p <= 1.6083e-19









% Compare individual correlation coefficients
p = compare_correlation_coefficients(r1,r2,n1,n2);

% Use cumulative of 1000 * N!







