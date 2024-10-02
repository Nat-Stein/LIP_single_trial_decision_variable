function [pValue,df,tStat] = pairedTforMeansSEM(data1,data2,SE1,SE2)
% perform paired t-test on pairs of \mu \pm s.e.
% Example paired data and standard errors
% 
% history. chatGPT3.5 provided this in response to Daniel Wolpert's
% conversation: three promtpts
% 1. is there a weighted t test
% 2. how can i do it in matlab [gave me an example where only one s.e per sample non paired]
% 3.can you give me a version where there are paired samples where  each sample has a standard error
% mns added the check for number of elements
%
% see testPairedTforMeansSEM.m 

if ~isscalar(unique([numel(data1), numel(data2), numel(SE1), numel(SE2)] ))
    error('all data and SE must contain the same number of elements')
end

% Calculate differences and combined standard error of differences
differences = data2 - data1;
combinedSE = sqrt(SE1.^2 + SE2.^2); % Assuming independence between pairs

% Calculate weights from the inverse of the variance of differences
weights = 1 ./ combinedSE.^2;

% Weighted mean of differences
weightedMeanDifferences = sum(weights .* differences) / sum(weights);

% Weighted variance of the differences
weightedVarDifferences = sum(weights .* (differences - weightedMeanDifferences).^2) / sum(weights);

% Standard error of the weighted mean difference
stdErrWeightedMeanDifference = sqrt(weightedVarDifferences / sum(weights));

% Perform a one-sample t-test on weighted differences
% Null hypothesis: mean difference = 0
tStat = weightedMeanDifferences / stdErrWeightedMeanDifference;
df = length(differences) - 1; % Degrees of freedom
pValue = 2 * tcdf(-abs(tStat), df); % Two-tailed test

% Display results
% fprintf('Weighted Mean of Differences: %f\n', weightedMeanDifferences);
% fprintf('t-Statistic: %f\n', tStat);
% fprintf('Degrees of Freedom: %f\n', df);
% fprintf('p-Value: %f\n', pValue);