function [avg_corr, sem_corr, z_scores, p, h] = averageCorrelation(correlation_list)
% Function to compute the average of correlation coefficients using Fisher's z transformation
% independently for each row
% Inputs:
%   - correlation_lists: A list of correlation coefficients
% Output:
%   - avg_corr: Average correlation coefficient after Fisher's z transformation
%   - sem_corr: Standard error of the mean of correlation coefficients

% Number of correlation coefficients
num_correlations = size(correlation_list,2);

% Convert correlation coefficients to Fisher's z scores
z_scores = 0.5 * log((1 + correlation_list) ./ (1 - correlation_list));

% Average the z scores
avg_z = sum(z_scores,2) / num_correlations;

% Convert the average z score back to correlation coefficient
avg_corr = (exp(2 * avg_z) - 1) ./ (exp(2 * avg_z) + 1);


% Compute the standard error of the mean
sem_corr = std(z_scores,[],2) / sqrt(num_correlations);

if size(correlation_list, 1)== 1
    [h, p] = ztest(0, avg_corr, sem_corr);
else
    h = 1; 
    p = NaN;
end

%[h_t,p_t,ci_t, tStat] = ttest(z_scores);


end
