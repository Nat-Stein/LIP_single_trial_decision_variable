

dimLabelsCS = {'ramp', 'PC1', 'Tin^{con}_{in}', 'What-decoder', 'When-decoder'};



%% Stats on cosine similarities
load(fullfile(saveLoc, 'cosine_sim'), 'cosine_sim')

allpVals = [];
for d1 = 1 : length(dimLabelsCS)
    for d2 = 1 : length(dimLabelsCS)
        if d1 < d2
            dVec = squeeze(cosine_sim(:, d1, d2));
            [h, p, tStat] = ttest(dVec);
            disp([dimLabelsCS{d1} ' vs ' dimLabelsCS{d2} ': p=' num2str(p)])
            disp(['mean +- sem = ' num2str(nanmean(dVec)) ' +- ' num2str(std(dVec)/sqrt(8))])
            allpVals = [allpVals p];
        end
    end
end
disp(['All p < ' num2str(max(allpVals))])

% Output
% ramp vs PC1: p=9.5073e-07
% mean +- sem = 0.64954 +- 0.040892
% ramp vs Tin: p=0.00045733
% mean +- sem = 0.34131 +- 0.055295
% ramp vs what decoder: p=0.00015175
% mean +- sem = 0.43956 +- 0.059545
% ramp vs when decoder: p=1.5045e-06
% mean +- sem = 0.62293 +- 0.041947
% PC1 vs Tin: p=0.00031195
% mean +- sem = 0.40499 +- 0.061611
% PC1 vs what decoder: p=0.00066464
% mean +- sem = 0.37741 +- 0.06509
% PC1 vs when decoder: p=0.00025631
% mean +- sem = 0.39929 +- 0.058835
% Tin vs what decoder: p=0.00022891
% mean +- sem = 0.28111 +- 0.040671
% Tin vs when decoder: p=4.2163e-05
% mean +- sem = 0.26063 +- 0.028907
% what decoder vs when decoder: p=0.0002512
% mean +- sem = 0.47409 +- 0.06963
% All p < 0.00066464

%% Stats on within-trial correlations

load(fullfile(saveLoc, 'resCorr_wiTrial'), 'filtW50b', 'resCorr', 'r_conv50b', 'rz_sess', 'all_zrMean');

allPs = [];
for dim1 = 1 : size(resCorr, 1)
    for dim2 = 1 : size(resCorr, 2)
        if dim1 == dim2
            rz_sess(dim1, dim2) = 1;
        else
            for sess =  1 : 8
                rvals = resCorr{dim1, dim2, sess};

                rvals_n = rvals(~isnan(rvals)); 
                [r_mean, sem_corr, z_scores, p] = averageCorrelation(rvals_n);
                % Same as:
                % dataOut = r2z(resCorr{dim1, dim2, sess, 1});
                % z_mean = nanmean(dataOut);
                % r_mean = tanh(z_mean);
                %all_zrMean{dim1, dim2}(sess) = r_mean; % This one
                [h,p,ci,stats] = ttest(z_scores);
                
                disp(['S' num2str(sess) ': ' dimLabelsCS{dim1} ' vs ' dimLabelsCS{dim2} ': p=' num2str(p)])
                allPs = [allPs p];
            end
            %rz_sess(dim1, dim2) = nanmean(all_zrMean{dim1, dim2}); % This one
        end
    end
end
disp(['All p <= ' num2str(max(allPs))])

% Output t-test on z-transpormed data
% S1: ramp vs PC1: p=0
% S2: ramp vs PC1: p=0
% S3: ramp vs PC1: p=0
% S4: ramp vs PC1: p=0
% S5: ramp vs PC1: p=0
% S6: ramp vs PC1: p=0
% S7: ramp vs PC1: p=0
% S8: ramp vs PC1: p=0
% S1: ramp vs Tin^{con}_{in}: p=0
% S2: ramp vs Tin^{con}_{in}: p=5.8328e-24
% S3: ramp vs Tin^{con}_{in}: p=0
% S4: ramp vs Tin^{con}_{in}: p=1.8709e-87
% S5: ramp vs Tin^{con}_{in}: p=0
% S6: ramp vs Tin^{con}_{in}: p=0
% S7: ramp vs Tin^{con}_{in}: p=1.3938e-231
% S8: ramp vs Tin^{con}_{in}: p=9.4199e-199
% S1: ramp vs What-decoder: p=9.6667e-296
% S2: ramp vs What-decoder: p=7.0757e-280
% S3: ramp vs What-decoder: p=0
% S4: ramp vs What-decoder: p=1.2632e-144
% S5: ramp vs What-decoder: p=2.4651e-144
% S6: ramp vs What-decoder: p=0
% S7: ramp vs What-decoder: p=1.5306e-297
% S8: ramp vs What-decoder: p=0
% S1: ramp vs When-decoder: p=0
% S2: ramp vs When-decoder: p=0
% S3: ramp vs When-decoder: p=0
% S4: ramp vs When-decoder: p=1.9643e-246
% S5: ramp vs When-decoder: p=0
% S6: ramp vs When-decoder: p=0
% S7: ramp vs When-decoder: p=0
% S8: ramp vs When-decoder: p=0
% S1: PC1 vs ramp: p=0
% S2: PC1 vs ramp: p=0
% S3: PC1 vs ramp: p=0
% S4: PC1 vs ramp: p=0
% S5: PC1 vs ramp: p=0
% S6: PC1 vs ramp: p=0
% S7: PC1 vs ramp: p=0
% S8: PC1 vs ramp: p=0
% S1: PC1 vs Tin^{con}_{in}: p=0
% S2: PC1 vs Tin^{con}_{in}: p=9.4057e-31
% S3: PC1 vs Tin^{con}_{in}: p=0
% S4: PC1 vs Tin^{con}_{in}: p=2.7348e-43
% S5: PC1 vs Tin^{con}_{in}: p=0
% S6: PC1 vs Tin^{con}_{in}: p=0
% S7: PC1 vs Tin^{con}_{in}: p=0
% S8: PC1 vs Tin^{con}_{in}: p=0
% S1: PC1 vs What-decoder: p=1.9763e-323
% S2: PC1 vs What-decoder: p=4.805e-105
% S3: PC1 vs What-decoder: p=0
% S4: PC1 vs What-decoder: p=7.1522e-57
% S5: PC1 vs What-decoder: p=1.6582e-109
% S6: PC1 vs What-decoder: p=0
% S7: PC1 vs What-decoder: p=0
% S8: PC1 vs What-decoder: p=0
% S1: PC1 vs When-decoder: p=0
% S2: PC1 vs When-decoder: p=1.5355e-184
% S3: PC1 vs When-decoder: p=0
% S4: PC1 vs When-decoder: p=1.2005e-62
% S5: PC1 vs When-decoder: p=0
% S6: PC1 vs When-decoder: p=0
% S7: PC1 vs When-decoder: p=0
% S8: PC1 vs When-decoder: p=0
% S1: Tin^{con}_{in} vs ramp: p=0
% S2: Tin^{con}_{in} vs ramp: p=5.8328e-24
% S3: Tin^{con}_{in} vs ramp: p=0
% S4: Tin^{con}_{in} vs ramp: p=1.8709e-87
% S5: Tin^{con}_{in} vs ramp: p=0
% S6: Tin^{con}_{in} vs ramp: p=0
% S7: Tin^{con}_{in} vs ramp: p=1.3938e-231
% S8: Tin^{con}_{in} vs ramp: p=9.4199e-199
% S1: Tin^{con}_{in} vs PC1: p=0
% S2: Tin^{con}_{in} vs PC1: p=9.4057e-31
% S3: Tin^{con}_{in} vs PC1: p=0
% S4: Tin^{con}_{in} vs PC1: p=2.7348e-43
% S5: Tin^{con}_{in} vs PC1: p=0
% S6: Tin^{con}_{in} vs PC1: p=0
% S7: Tin^{con}_{in} vs PC1: p=0
% S8: Tin^{con}_{in} vs PC1: p=0
% S1: Tin^{con}_{in} vs What-decoder: p=0
% S2: Tin^{con}_{in} vs What-decoder: p=3.1915e-89
% S3: Tin^{con}_{in} vs What-decoder: p=0
% S4: Tin^{con}_{in} vs What-decoder: p=1.2291e-67
% S5: Tin^{con}_{in} vs What-decoder: p=2.798e-106
% S6: Tin^{con}_{in} vs What-decoder: p=0
% S7: Tin^{con}_{in} vs What-decoder: p=0
% S8: Tin^{con}_{in} vs What-decoder: p=2.0442e-301
% S1: Tin^{con}_{in} vs When-decoder: p=0
% S2: Tin^{con}_{in} vs When-decoder: p=3.112e-74
% S3: Tin^{con}_{in} vs When-decoder: p=0
% S4: Tin^{con}_{in} vs When-decoder: p=1.7316e-46
% S5: Tin^{con}_{in} vs When-decoder: p=0
% S6: Tin^{con}_{in} vs When-decoder: p=0
% S7: Tin^{con}_{in} vs When-decoder: p=0
% S8: Tin^{con}_{in} vs When-decoder: p=4.5894e-310
% S1: What-decoder vs ramp: p=9.6667e-296
% S2: What-decoder vs ramp: p=7.0757e-280
% S3: What-decoder vs ramp: p=0
% S4: What-decoder vs ramp: p=1.2632e-144
% S5: What-decoder vs ramp: p=2.4651e-144
% S6: What-decoder vs ramp: p=0
% S7: What-decoder vs ramp: p=1.5306e-297
% S8: What-decoder vs ramp: p=0
% S1: What-decoder vs PC1: p=1.9763e-323
% S2: What-decoder vs PC1: p=4.805e-105
% S3: What-decoder vs PC1: p=0
% S4: What-decoder vs PC1: p=7.1522e-57
% S5: What-decoder vs PC1: p=1.6582e-109
% S6: What-decoder vs PC1: p=0
% S7: What-decoder vs PC1: p=0
% S8: What-decoder vs PC1: p=0
% S1: What-decoder vs Tin^{con}_{in}: p=0
% S2: What-decoder vs Tin^{con}_{in}: p=3.1915e-89
% S3: What-decoder vs Tin^{con}_{in}: p=0
% S4: What-decoder vs Tin^{con}_{in}: p=1.2291e-67
% S5: What-decoder vs Tin^{con}_{in}: p=2.798e-106
% S6: What-decoder vs Tin^{con}_{in}: p=0
% S7: What-decoder vs Tin^{con}_{in}: p=0
% S8: What-decoder vs Tin^{con}_{in}: p=2.0442e-301
% S1: What-decoder vs When-decoder: p=0
% S2: What-decoder vs When-decoder: p=0
% S3: What-decoder vs When-decoder: p=0
% S4: What-decoder vs When-decoder: p=1.0945e-145
% S5: What-decoder vs When-decoder: p=4.2388e-244
% S6: What-decoder vs When-decoder: p=0
% S7: What-decoder vs When-decoder: p=0
% S8: What-decoder vs When-decoder: p=0
% S1: When-decoder vs ramp: p=0
% S2: When-decoder vs ramp: p=0
% S3: When-decoder vs ramp: p=0
% S4: When-decoder vs ramp: p=1.9643e-246
% S5: When-decoder vs ramp: p=0
% S6: When-decoder vs ramp: p=0
% S7: When-decoder vs ramp: p=0
% S8: When-decoder vs ramp: p=0
% S1: When-decoder vs PC1: p=0
% S2: When-decoder vs PC1: p=1.5355e-184
% S3: When-decoder vs PC1: p=0
% S4: When-decoder vs PC1: p=1.2005e-62
% S5: When-decoder vs PC1: p=0
% S6: When-decoder vs PC1: p=0
% S7: When-decoder vs PC1: p=0
% S8: When-decoder vs PC1: p=0
% S1: When-decoder vs Tin^{con}_{in}: p=0
% S2: When-decoder vs Tin^{con}_{in}: p=3.112e-74
% S3: When-decoder vs Tin^{con}_{in}: p=0
% S4: When-decoder vs Tin^{con}_{in}: p=1.7316e-46
% S5: When-decoder vs Tin^{con}_{in}: p=0
% S6: When-decoder vs Tin^{con}_{in}: p=0
% S7: When-decoder vs Tin^{con}_{in}: p=0
% S8: When-decoder vs Tin^{con}_{in}: p=4.5894e-310
% S1: When-decoder vs What-decoder: p=0
% S2: When-decoder vs What-decoder: p=0
% S3: When-decoder vs What-decoder: p=0
% S4: When-decoder vs What-decoder: p=1.0945e-145
% S5: When-decoder vs What-decoder: p=4.2388e-244
% S6: When-decoder vs What-decoder: p=0
% S7: When-decoder vs What-decoder: p=0
% S8: When-decoder vs What-decoder: p=0
% All p <= 5.8328e-24



% Output ztest
% S1: ramp vs PC1: p=0
% S2: ramp vs PC1: p=0
% S3: ramp vs PC1: p=0
% S4: ramp vs PC1: p=0
% S5: ramp vs PC1: p=0
% S6: ramp vs PC1: p=0
% S7: ramp vs PC1: p=0
% S8: ramp vs PC1: p=0
% S1: ramp vs Tin^{con}_{in}: p=0
% S2: ramp vs Tin^{con}_{in}: p=1.6239e-24
% S3: ramp vs Tin^{con}_{in}: p=0
% S4: ramp vs Tin^{con}_{in}: p=3.173e-94
% S5: ramp vs Tin^{con}_{in}: p=0
% S6: ramp vs Tin^{con}_{in}: p=0
% S7: ramp vs Tin^{con}_{in}: p=2.8252e-260
% S8: ramp vs Tin^{con}_{in}: p=1.4505e-218
% S1: ramp vs What-decoder: p=0
% S2: ramp vs What-decoder: p=0
% S3: ramp vs What-decoder: p=0
% S4: ramp vs What-decoder: p=3.8476e-165
% S5: ramp vs What-decoder: p=1.6152e-160
% S6: ramp vs What-decoder: p=0
% S7: ramp vs What-decoder: p=0
% S8: ramp vs What-decoder: p=0
% S1: ramp vs When-decoder: p=0
% S2: ramp vs When-decoder: p=0
% S3: ramp vs When-decoder: p=0
% S4: ramp vs When-decoder: p=1.1702e-310
% S5: ramp vs When-decoder: p=0
% S6: ramp vs When-decoder: p=0
% S7: ramp vs When-decoder: p=0
% S8: ramp vs When-decoder: p=0
% S1: PC1 vs ramp: p=0
% S2: PC1 vs ramp: p=0
% S3: PC1 vs ramp: p=0
% S4: PC1 vs ramp: p=0
% S5: PC1 vs ramp: p=0
% S6: PC1 vs ramp: p=0
% S7: PC1 vs ramp: p=0
% S8: PC1 vs ramp: p=0
% S1: PC1 vs Tin^{con}_{in}: p=0
% S2: PC1 vs Tin^{con}_{in}: p=1.0946e-31
% S3: PC1 vs Tin^{con}_{in}: p=0
% S4: PC1 vs Tin^{con}_{in}: p=9.7625e-45
% S5: PC1 vs Tin^{con}_{in}: p=0
% S6: PC1 vs Tin^{con}_{in}: p=0
% S7: PC1 vs Tin^{con}_{in}: p=0
% S8: PC1 vs Tin^{con}_{in}: p=0
% S1: PC1 vs What-decoder: p=0
% S2: PC1 vs What-decoder: p=2.5517e-116
% S3: PC1 vs What-decoder: p=0
% S4: PC1 vs What-decoder: p=1.3908e-59
% S5: PC1 vs What-decoder: p=1.8133e-119
% S6: PC1 vs What-decoder: p=0
% S7: PC1 vs What-decoder: p=0
% S8: PC1 vs What-decoder: p=0
% S1: PC1 vs When-decoder: p=0
% S2: PC1 vs When-decoder: p=1.5583e-215
% S3: PC1 vs When-decoder: p=0
% S4: PC1 vs When-decoder: p=2.2759e-66
% S5: PC1 vs When-decoder: p=0
% S6: PC1 vs When-decoder: p=0
% S7: PC1 vs When-decoder: p=0
% S8: PC1 vs When-decoder: p=0
% S1: Tin^{con}_{in} vs ramp: p=0
% S2: Tin^{con}_{in} vs ramp: p=1.6239e-24
% S3: Tin^{con}_{in} vs ramp: p=0
% S4: Tin^{con}_{in} vs ramp: p=3.173e-94
% S5: Tin^{con}_{in} vs ramp: p=0
% S6: Tin^{con}_{in} vs ramp: p=0
% S7: Tin^{con}_{in} vs ramp: p=2.8252e-260
% S8: Tin^{con}_{in} vs ramp: p=1.4505e-218
% S1: Tin^{con}_{in} vs PC1: p=0
% S2: Tin^{con}_{in} vs PC1: p=1.0946e-31
% S3: Tin^{con}_{in} vs PC1: p=0
% S4: Tin^{con}_{in} vs PC1: p=9.7625e-45
% S5: Tin^{con}_{in} vs PC1: p=0
% S6: Tin^{con}_{in} vs PC1: p=0
% S7: Tin^{con}_{in} vs PC1: p=0
% S8: Tin^{con}_{in} vs PC1: p=0
% S1: Tin^{con}_{in} vs What-decoder: p=0
% S2: Tin^{con}_{in} vs What-decoder: p=4.3373e-97
% S3: Tin^{con}_{in} vs What-decoder: p=0
% S4: Tin^{con}_{in} vs What-decoder: p=5.4873e-72
% S5: Tin^{con}_{in} vs What-decoder: p=2.1454e-116
% S6: Tin^{con}_{in} vs What-decoder: p=0
% S7: Tin^{con}_{in} vs What-decoder: p=0
% S8: Tin^{con}_{in} vs What-decoder: p=0
% S1: Tin^{con}_{in} vs When-decoder: p=0
% S2: Tin^{con}_{in} vs When-decoder: p=7.3571e-80
% S3: Tin^{con}_{in} vs When-decoder: p=0
% S4: Tin^{con}_{in} vs When-decoder: p=1.8473e-48
% S5: Tin^{con}_{in} vs When-decoder: p=0
% S6: Tin^{con}_{in} vs When-decoder: p=0
% S7: Tin^{con}_{in} vs When-decoder: p=0
% S8: Tin^{con}_{in} vs When-decoder: p=0
% S1: What-decoder vs ramp: p=0
% S2: What-decoder vs ramp: p=0
% S3: What-decoder vs ramp: p=0
% S4: What-decoder vs ramp: p=3.8476e-165
% S5: What-decoder vs ramp: p=1.6152e-160
% S6: What-decoder vs ramp: p=0
% S7: What-decoder vs ramp: p=0
% S8: What-decoder vs ramp: p=0
% S1: What-decoder vs PC1: p=0
% S2: What-decoder vs PC1: p=2.5517e-116
% S3: What-decoder vs PC1: p=0
% S4: What-decoder vs PC1: p=1.3908e-59
% S5: What-decoder vs PC1: p=1.8133e-119
% S6: What-decoder vs PC1: p=0
% S7: What-decoder vs PC1: p=0
% S8: What-decoder vs PC1: p=0
% S1: What-decoder vs Tin^{con}_{in}: p=0
% S2: What-decoder vs Tin^{con}_{in}: p=4.3373e-97
% S3: What-decoder vs Tin^{con}_{in}: p=0
% S4: What-decoder vs Tin^{con}_{in}: p=5.4873e-72
% S5: What-decoder vs Tin^{con}_{in}: p=2.1454e-116
% S6: What-decoder vs Tin^{con}_{in}: p=0
% S7: What-decoder vs Tin^{con}_{in}: p=0
% S8: What-decoder vs Tin^{con}_{in}: p=0
% S1: What-decoder vs When-decoder: p=0
% S2: What-decoder vs When-decoder: p=0
% S3: What-decoder vs When-decoder: p=0
% S4: What-decoder vs When-decoder: p=1.809e-167
% S5: What-decoder vs When-decoder: p=1.0073e-278
% S6: What-decoder vs When-decoder: p=0
% S7: What-decoder vs When-decoder: p=0
% S8: What-decoder vs When-decoder: p=0
% S1: When-decoder vs ramp: p=0
% S2: When-decoder vs ramp: p=0
% S3: When-decoder vs ramp: p=0
% S4: When-decoder vs ramp: p=1.1702e-310
% S5: When-decoder vs ramp: p=0
% S6: When-decoder vs ramp: p=0
% S7: When-decoder vs ramp: p=0
% S8: When-decoder vs ramp: p=0
% S1: When-decoder vs PC1: p=0
% S2: When-decoder vs PC1: p=1.5583e-215
% S3: When-decoder vs PC1: p=0
% S4: When-decoder vs PC1: p=2.2759e-66
% S5: When-decoder vs PC1: p=0
% S6: When-decoder vs PC1: p=0
% S7: When-decoder vs PC1: p=0
% S8: When-decoder vs PC1: p=0
% S1: When-decoder vs Tin^{con}_{in}: p=0
% S2: When-decoder vs Tin^{con}_{in}: p=7.3571e-80
% S3: When-decoder vs Tin^{con}_{in}: p=0
% S4: When-decoder vs Tin^{con}_{in}: p=1.8473e-48
% S5: When-decoder vs Tin^{con}_{in}: p=0
% S6: When-decoder vs Tin^{con}_{in}: p=0
% S7: When-decoder vs Tin^{con}_{in}: p=0
% S8: When-decoder vs Tin^{con}_{in}: p=0
% S1: When-decoder vs What-decoder: p=0
% S2: When-decoder vs What-decoder: p=0
% S3: When-decoder vs What-decoder: p=0
% S4: When-decoder vs What-decoder: p=1.809e-167
% S5: When-decoder vs What-decoder: p=1.0073e-278
% S6: When-decoder vs What-decoder: p=0
% S7: When-decoder vs What-decoder: p=0
% S8: When-decoder vs What-decoder: p=0
% All p <= 1.6239e-24








