% calc_wiTrialCorr
% Used in Figure 4e

clear r_conv50b
for sess = 1 : 8
    isSess = find(sigSL.session == sess); % All trials in this session
    sig_coh_sess = sigSL.sig_coh(isSess); % signed coherences of all trials in this session
    cohs = unique(sig_coh_sess);
    for d = 1 : length(useD)
        dim = dNames{useD(d)}; % name of direction in state space
        assignSignal % data in direction in state space (dim)
        % Create residuals
        dMat = sl(isSess, :);
        if size(dMat, 1) > 0
            rMat = nan(size(dMat));
            for c = 1 : length(cohs)
                coh = cohs(c);
                trls = find(sig_coh_sess == coh);
                
                rMat(trls, :) = dMat(trls, :) - repmat(nanmean(dMat(trls, :), 1), length(trls), 1); % Compute residual within coherence
            end
            r_conv50b{d, sess} = conv2(1, filtW50b, rMat, 'same'); % filtW50b: 50ms boxcar filter defined in general_settings
        end
    end
end

rtAdd = 0.1;
for sess = 1 : 8
    isSess = find(sigSL.session == sess); % All trials in this session
    rt = sigSL.rt(isSess);
    for dim1 = 1 : size(r_conv50b, 1)
        for dim2 = 1 : size(r_conv50b, 1)
            % par = par_sess{sess};
            if size(r_conv50b{dim2, sess}, 1) > 0 && size(r_conv50b{dim1, sess}, 1) > 0
                resCorr{dim1, dim2, sess} = nan(1, size(r_conv50b{dim2, sess}, 1));
                for trl = 1 : size(r_conv50b{dim2, sess}, 1)
                    tpts = find(~isnan(r_conv50b{dim1, sess}(trl, :)) & ~isnan(r_conv50b{dim2, sess}(trl, :)) &...
                        sigSL.t >= 0.2 & sigSL.t <= rt(trl) - rtAdd);
                    
                    if length(tpts) > 5
                        tpts2 = tpts(1 : 50 : end);
                        if length(tpts2) >= 4
                            [R, P] = corrcoef(r_conv50b{dim1, sess}(trl, tpts2), r_conv50b{dim2, sess}(trl, tpts2), 'rows', 'complete');
                            resCorr{dim1, dim2, sess}(trl) = R(1, 2);
                        end
                    end
                    
                end
            end
        end
    end
end


%% Extract the most important values and sort by dimesnion order
clear allMeds all_zrMean rz_sess
for dim1 = 1 : size(resCorr, 1)
    for dim2 = 1 : size(resCorr, 2)
        if dim1 == dim2
            rz_sess(dim1, dim2) = 1;
        else
            for sess =  1 : 8
                d1 = resCorr{dim1, dim2, sess};

                d1_n = d1(~isnan(d1)); r_mean = averageCorrelation(d1_n);  % This one
                % Same as:
                % dataOut = r2z(resCorr{dim1, dim2, sess, 1});
                % z_mean = nanmean(dataOut);
                % r_mean = tanh(z_mean);
                all_zrMean{dim1, dim2}(sess) = r_mean; % This one
            end
            rz_sess(dim1, dim2) = nanmean(all_zrMean{dim1, dim2}); % This one
        end
    end
end

save(fullfile(saveLoc, 'resCorr_wiTrial'), 'par_sess', 'filtW50b', 'resCorr', 'r_conv50b', 'rz_sess', 'all_zrMean', '-v7.3');



