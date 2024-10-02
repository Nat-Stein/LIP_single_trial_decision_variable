%% Heatmaps of cosine similarity & within-trial correlation of diffusion across signals

% -------------------------------------------------------
% Calculate within-trial correlation
fw = 80; x = 0:1:fw; y = gaussmf(x,[fw/4 fw/2]);
filtW = y / sum(y);
fw = 50; x = 0:1:fw; y = gaussmf(x,[fw/4 fw/2]);
filtW50 = y / sum(y);
fw = 50; filtW50b = ones(1, fw)/fw;

% mainDim = [1]; % bioRxiv 
% compDims = [3 26 20 24 31 32 57 58 59 60 61]; % bioRxiv 
mainDim = [63 20 1 61 7 62]; % change to match input of sigSL
compDims = mainDim; 

doSave = 1; 

calc_wiTrialCorr % old version: calc_withinTrialCorr

% -------------------------------------------------------
% Calculate cosine similarity

calc_weightCS

% save('E:\Matalb analyses\cosine_sim', 'cosine_sim')
save('E:\Matalb analyses\allCs', 'allCs')

% -------------------------------------------------------
% Plot matrices
% load('E:\Matalb analyses\resCorr_wiTrial_230427')
load('E:\Matalb analyses\resCorr_wiTrial_230522')
fullfile(saveLoc, 'resCorr_wiTrial')
load('E:\Matalb analyses\darkRedBlue')
load('E:\Matalb analyses\allCs')

mainDim = [63 20 1 61 7];
compDims = mainDim;
dimLabels = {'ramp', 'PC1', 'Tin', 'what decoder', 'when decoder'};
% dimLabels = {'ramp', 'ramp no lasso'};

clear allMeds all_zrMean
for dim1 = mainDim
    for dim2 = compDims
        if dim1 == dim2
            r_sess(dim1, dim2) = 1;
            rz_sess(dim1, dim2) = 1;
        else
            for sess =  1 : 8
                d1 = resCorr_wiTrialConv50bStep{dim1, dim2, sess, 1}; % Filtered

                d1_n = d1(~isnan(d1)); [r_mean, sem_corr] = averageCorrelation(d1_n);
                % Same as:
                % dataOut = r2z(resCorr_wiTrialConv50bStep{dim1, dim2, sess, 1});
                % z_mean = nanmean(dataOut);
                % r_mean = tanh(z_mean);
                all_zrMean{dim1, dim2}(sess) = r_mean;
                
                %                 d1 = resCorr_wiTrial{dim1, dim2, sess, 1};
                % d1 = meanCorrZ{dim1, dim2, sess, 1};
                useTrl = find(~isnan(d1));
                med1 = median(d1(useTrl));
                allMeds{dim1, dim2}(sess) = med1;
                
%                 allMeds{dim1, dim2}(sess) = meanCorrR{dim1, dim2, sess, 1};
            end
            r_sess(dim1, dim2) = nanmean(allMeds{dim1, dim2});
            rz_sess(dim1, dim2) = nanmean(all_zrMean{dim1, dim2}); % This one
%             r_sessMean(dim1, dim2) = nanmean(allMeds{dim1, dim2});
        end
    end
end
% Correlation between Ramp and RampLasso: rz_sess = 0.9912


%% Plotting
% Turn into: plot_Fig4d % not written yet

useDims = flip(mainDim);
load('E:\Matalb analyses\darkRedBlue')
figure('Position', [600 400 800 400]); hold on;
subplot(1, 2, 1); hold on; title('r')
A = rz_sess(useDims, useDims); % A = r_sess(useDims, useDims); 
A = triu(A);
imagesc(A, [-1 1]);
set(gca,'xtick',1:length(dimLabels),'xticklabel', flip(dimLabels),'tickdir','out');
set(gca,'ytick',1:length(dimLabels),'yticklabel', flip(dimLabels),'tickdir','out');
xtickangle(45)
colormap(darkRedBlue);
colorbar
axis tight
axis square


subplot(1, 2, 2); hold on; title('CS')
A = squeeze(nanmean(allCs(1:8, useDims, useDims), 1)); A = triu(A);
imagesc(A, [-1 1]);
set(gca,'xtick',1:length(dimLabels),'xticklabel', flip(dimLabels),'tickdir','out');
set(gca,'ytick',1:length(dimLabels),'yticklabel', flip(dimLabels),'tickdir','out');
xtickangle(45)
colormap(darkRedBlue);
colorbar
axis tight
axis square

%% -------------------------------------------------------
% Stats
% Do stats on CS
allVals = [];
d1 = 0;
for dim1 = mainDim
    d1 = d1+1; d2=0;
    for dim2 = compDims
        d2 = d2 +1;
        if dim1 < dim2
            dVec = allCs(1:8, dim1, dim2);
            [t, p] = ttest(dVec);
            disp([dimLabels{d1} ' vs ' dimLabels{d2} ': t(7)=' num2str(t) ', p=' num2str(p)])
            disp(['mean +- sem = ' num2str(nanmean(dVec)) ' +- ' num2str(std(dVec)/sqrt(8))])
            allVals = [allVals squeeze(allCs(1:8, dim1, dim2))'];
        end
    end
end
[t, p] = ttest(allVals);
disp(['all CS: t(' num2str(length(allVals)-1) ')=' num2str(t) ', p=' num2str(p)])

% CS stats
% PC1 vs Ramp: t(7)=1, p=3.1945e-07
% PC1 vs what-decoder: t(7)=1, p=0.002794
% Tin_new vs Ramp: t(7)=1, p=0.00023187
% Tin_new vs PC1: t(7)=1, p=0.00014517
% Tin_new vs what-decoder: t(7)=1, p=0.00023097
% Tin_new vs reg2whenC: t(7)=1, p=4.2163e-05
% what-decoder vs Ramp: t(7)=1, p=0.0032425
% reg2whenC vs Ramp: t(7)=1, p=3.0887e-08
% reg2whenC vs PC1: t(7)=1, p=0.0002917
% reg2whenC vs what-decoder: t(7)=1, p=0.00082971
% all CS: t(79)=1, p=2.0541e-29
% Old
% PC1 vs Ramp: t(7)=1, p=4.5229e-07
% PC1 vs what-decoder: t(7)=1, p=0.002794
% Tin_new vs Ramp: t(7)=1, p=0.00013799
% Tin_new vs PC1: t(7)=1, p=0.00014517
% Tin_new vs what-decoder: t(7)=1, p=0.00023097
% Tin_new vs reg2whenC: t(7)=1, p=4.2163e-05
% what-decoder vs Ramp: t(7)=1, p=0.0048116
% reg2whenC vs Ramp: t(7)=1, p=2.2243e-08
% reg2whenC vs PC1: t(7)=1, p=0.0002917
% reg2whenC vs what-decoder: t(7)=1, p=0.00082971
% all CS: t(79)=1, p=2.0215e-29


nanmean(allCs(:, 7, 61))
nanstd(allCs(:, 7, 61))


% -------------------------------------------------------
% Do stats on correlations? 
[avg_corr, sem_corr, z_scores] = averageCorrelation(correlation_list);

allVals = [];
for dim1 = mainDim
    for dim2 = compDims
        if dim1 < dim2
            [t, p] = ttest(allCs(1:8, dim1, dim2));
            disp([dimName{dim1} ' vs ' dimName{dim2} ': t(7)=' num2str(t) ', p=' num2str(p)])
            allVals = [allVals squeeze(allCs(1:8, dim1, dim2))'];
        end
    end
end
[t, p] = ttest(allVals)
disp(['all CS: t(' num2str(length(allVals)-1) ')=' num2str(t) ', p=' num2str(p)])
















