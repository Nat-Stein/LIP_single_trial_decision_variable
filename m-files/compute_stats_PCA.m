% compute_stats_PCA

inclPC2 = 100;
allExpl = nan(8, inclPC2);
for sess = 1 : 8
    clear outPCA
    suff = '_mean_stp50'; load(fullfile(saveLoc, ['pca' suff '_S' num2str(sess)]), 'outPCA');
    numPCs(sess) = length(outPCA.explainedD); 
    nn = min(numPCs(sess), inclPC2);
    allExpl(sess, 1 : nn) = outPCA.explainedD(1:nn);
    partRat(sess) = sum(outPCA.latentD)^2 / sum(outPCA.latentD.^2);
end

allPC_Expl = cumsum(allExpl');
meanPC_Expl = mean(allPC_Expl, 2); 
% expl90 = find(meanPC_Expl >= 90); expl90 = expl90(1);
for sess = 1 : 8
    ex90 = find(allPC_Expl(:, sess) >= 90); ex90 = ex90(1);
    expl90_sess(sess) = ex90;
end
figure('Position', [200 200 300 300]); hold on
expl90 = find(meanPC_Expl >= 90); expl90 = expl90(1);
plot(allPC_Expl, 'Color', 0.2 * [1 1 1]); plot(meanPC_Expl, 'k', 'LineWidth', 2); 
plot(expl90, meanPC_Expl(expl90), 'ro', 'MarkerFaceColor','r')
plot(3, meanPC_Expl(3), 'bo', 'MarkerFaceColor','b')
ylim([0 100])
xlabel('PCs included')
ylabel('Cumulative variance explained')
set(gca,'TickDir','out', 'FontSize', fontSize);
disp(['Participation ratio (mean +- sem) = ' num2str(round(nanmean(partRat), 1))...
    ' +- ' num2str(round(nanstd(partRat)/sqrt(8), 1))])
disp(['90% of variance explained by (mean +- sem) = ' num2str(round(nanmean(expl90_sess), 1))...
    ' +- ' num2str(round(nanstd(expl90_sess)/sqrt(8), 1)) ' components'])
disp(['First PC explains ' num2str(round(nanmean(allPC_Expl(1, :)), 1))...
    ' +- ' num2str(round(nanstd(allPC_Expl(1, :))/sqrt(8), 1)) ' % of the variance  (mean +- std)'])
disp(['First 3 PCs explain ' num2str(round(nanmean(allPC_Expl(3, :)), 1))...
    ' +- ' num2str(round(nanstd(allPC_Expl(3, :))/sqrt(8), 1)) ' % of the variance  (mean +- std)'])

% Participation ratio (mean +- sem) = 4.5 +- 0.4
% 90% of variance explained by (mean +- sem) = 22.3 +- 2.6 components
% First PC explains 44 +- 2.5 % of the variance  (mean +- std)
% First 3 PCs explain 67.1 +- 3.1 % of the variance  (mean +- std)



