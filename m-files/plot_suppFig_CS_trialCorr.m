
% Load data

load(fullfile(saveLoc, 'cosine_sim'), 'cosine_sim')
load(fullfile(saveLoc, 'resCorr_wiTrial'), 'all_zrMean');
load(fullfile(saveLoc, 'cosine_sim_randperm'), 'cosine_sim_randperm') % cosine_sim_randperm(sess, d1, d2, :)
load(fullfile(saveLoc, 'resCorr_wiTrial_rand_gauss'), 'par_sess', 'filtW50b', 'rz_sess_rand') % , 'cosine_sim_rand');
load(fullfile(saveLoc, 'cosine_sim_rand'), 'cosine_sim_rand')
% load(fullfile(saveLoc, 'resCorr_rand_permute'), 'resCorr_rand_permute')
load(fullfile(saveLoc, 'single_trial_corr'))

dimLabelsCS = {'ramp', 'PC1', 'T^{con}_{in}', 'What-decoder', 'When-decoder'};



plot_rand = 1;
plot_perm = 1; 
% load('E:\Matalb analyses\darkRedBlue')
useD = [3 4 1 5 6];
colms = length(useD); 
if plot_rand == 1; colms = colms + 1; end
rows = 2;
mSize = 5;
compCol = [173 130 12]/255;
compCol2 = [54 38 191]/255;
ylims = [-0.1 1];

figure('Position', [600 400 800 400]); hold on;
for d = 1 : length(useD)
    dim1 = (d);
    pc = 0;
    legs = [];
    for d2 = 1 : length(useD)
        dim2 = (d2);
        if d ~= d2
            pc = pc + 1;
            subplot(rows, colms, d); hold on; % title('Cosine similarity');
            title(dimLabelsCS{d})
            % Real data
            cs = squeeze(cosine_sim(:, dim1, dim2));
            m = nanmean(cs);
            sem = nanstd(cs)/sqrt(sum(~isnan(cs)));
            plot(pc, m, 'ko', 'MarkerFaceColor','k', 'MarkerSize', mSize)
            line(pc*[1 1], [m-sem m+sem], 'Color', 'k')
            if plot_perm == 1
                % Permuted weights
                thisCol = compCol;
                cs = squeeze(nanmean(cosine_sim_randperm(:, dim1, dim2, :), 1));
                m = nanmean(cs);
                sem = nanstd(cs);
                plot(pc, m, 'o', 'Color', thisCol, 'MarkerFaceColor',thisCol, 'MarkerSize', mSize)
                line(pc*[1 1], [m-sem m+sem], 'Color', thisCol)
            end
            ylim(ylims)
            xlim([0 5])
            if d == 1; ylabel('cosine similarity'); end
            
            subplot(rows, colms, d + colms);  hold on; % title('r')
            
            cs = all_zrMean{dim1, dim2};
            m = nanmean(cs);
            sem = nanstd(cs)/sqrt(sum(~isnan(cs)));
            plot(pc, m, 'ko', 'MarkerFaceColor','k', 'MarkerSize', mSize)
            line(pc*[1 1], [m-sem m+sem], 'Color', 'k')
            if plot_perm == 1
                % Permuted weights
                thisCol = compCol;
                % cs = squeeze(nanmean(cosine_sim_randperm(:, dim1, dim2, :), 1));
                if dim2 > dim1
                    cs = nanmean(resCorr_rand_permute{:, dim2, dim1}, 2);
                else
                    cs = nanmean(resCorr_rand_permute{:, dim1, dim2}, 2);
                end
                m = nanmean(cs);
                sem = nanstd(cs);
                plot(pc, m, 'o', 'Color', thisCol, 'MarkerFaceColor',thisCol, 'MarkerSize', mSize)
                line(pc*[1 1], [m-sem m+sem], 'Color', thisCol)
            end
            ylim(ylims)
            xlim([0 5])
            if d == 1; ylabel('within-trial correlation'); end
            legs{pc} = dimLabelsCS{d2};
        end
    end
    set(gca,'xtick',1:length(legs),'xticklabel', legs,'tickdir','out');
    xtickangle(45)
    set(gca,'TickDir','out', 'FontSize', fontSize);
    subplot(rows, colms, d); hold on; 
    set(gca,'xtick',1:length(legs),'xticklabel', {'','','',''},'tickdir','out');
    set(gca,'TickDir','out', 'FontSize', fontSize);
end


if plot_rand == 1
    pc = 1; 
    subplot(rows, colms, d+1); hold on;
    title('random')
    % Random vectors - CS
    thisCol = compCol2;
    cs = squeeze(nanmean(cosine_sim_rand, 1));
    m = nanmean(cs);
    sem = nanstd(cs);
    plot(pc, m, 'o', 'Color', thisCol, 'MarkerFaceColor',thisCol, 'MarkerSize', mSize)
    line(pc*[1 1], [m-sem m+sem], 'Color', thisCol)
    ylim(ylims)
    set(gca,'xtick',pc,'xticklabel', {''},'tickdir','out');
    set(gca,'TickDir','out', 'FontSize', fontSize);
    
    subplot(rows, colms, d + 1 + colms);  hold on;
    % Random vectors - within-trial correlation
    thisCol = compCol2;
    cs = squeeze(rz_sess_rand);
    m = nanmean(cs);
    sem = nanstd(cs);
    plot(pc, m, 'o', 'Color', thisCol, 'MarkerFaceColor',thisCol, 'MarkerSize', mSize)
    line(pc*[1 1], [m-sem m+sem], 'Color', thisCol)
    ylim(ylims)
    set(gca,'xtick',pc ,'xticklabel', 'random','tickdir','out');
    xtickangle(45)
    set(gca,'TickDir','out', 'FontSize', fontSize);
    
end

