% Compute the When-decoder

%% Adaptend from AZ 2022
global par_sess
ndatasets = length(par_sess);

plot_figs = 0;

names = {'all','Tin','Tout','Min_C','Min_i','allButTin','allButTinOut'};
AUC = nan(ndatasets, length(names));

clear allWhenDs
for sess = 1 : ndatasets
    
    disp(num2str(sess))
    
    par = par_sess{sess};
    
    dt_ms = 25;
    dataset = sess;
    disp('Loading...')
    clear dat
    dat = get_data(dt_ms, [], dataset, par, 'dt_rel_RT', 0.05);
    
    struct2vars(dat);
    
    idx_include = choice==0;
    % idx_include = choice==1;
    
    dt = t(2)-t(1);
    
    R = [];
    D = [];
    T = [];
    trial_id = [];
    for i=1:length(t)-1
        r = zeros(size(RT));
        r(RT>=t(i) & RT<(t(i) + 0.15)) = 1;
        R = [R; r(idx_include)];
        T = [T; t(i)*ones(sum(idx_include),1)];
        
        D = cat(1,D,squeeze(H(:,i,idx_include))');
        
        trial_id = [trial_id; find(idx_include)];
    end
    
    % remove nans
    I = all(isnan(D),2);
    R = R(~I);
    D = D(~I,:);
    T = T(~I,:);
    trial_id = trial_id(~I);
    
    % sort by trial, easier to plot
    [~,ind] = sort(trial_id);
    trial_id = trial_id(ind);
    R = R(ind);
    D = D(ind,:);
    T = T(ind,:);
    
    %% regress
    
    allneurons = [1:dat.nneurons]';
    idx_neurons = {allneurons, dat.unitIdxLIP_Tin, dat.unitIdxLIP_Tout, dat.minCells.DinRFc, dat.minCells.DinRFi,...
        find(~ismember(allneurons, dat.unitIdxLIP_Tin)), find(~ismember(allneurons, [dat.unitIdxLIP_Tin, dat.unitIdxLIP_Tout]))};
    
    J = iseven(trial_id);
    W = zeros(dat.nneurons,length(idx_neurons));
    for k=1:length(idx_neurons)
        
        if ~isempty(idx_neurons{k})
            
            [B,FitInfo] = lassoglm(D(J,idx_neurons{k}),R(J),'binomial','link','logit','lambda',0.01);
            B0 = FitInfo.Intercept;
            
            coef = [B0; B];
            
            % pred
            yhat = glmval(coef,D(:,idx_neurons{k}),'logit');
            
            %% calc ROC
            scores = yhat(~J);
            labels = R(~J);
            posclass = 1;
            [~,~,~,auc] = perfcurve(labels,scores,posclass);
            AUC(sess,k) = auc;
            
            %%
            W(idx_neurons{k},k) = B;
            
            
            %% Plot the fit for Session 6 with all neurons (Supp Fig 6)
            if k==1 && sess == 6
                
                ntr_for_plot = 8;
                first_trial = 350;
                
                tr = unique(trial_id);
                tr_ind = findclose(tr, first_trial);
                tr_ind = tr(tr_ind:tr_ind + ntr_for_plot - 1);
                
                figure; hold on
                set(gcf,'Position',[432  572  707  313]);
                colores = movshon_colors(3);
                
                vR = [];
                vYhat = [];
                for i=1:length(tr_ind)
                    ind = trial_id ==tr_ind(i) & T>0.2;
                    vR = [vR; R(ind); NaN];
                    vYhat = [vYhat; yhat(ind); NaN];
                end
                
                tt = dt * [1:length(vR)] - dt;
                hs = stairs(tt, vR,'k');
                set(hs,'LineWidth',1);
                hold all
                plot(tt, vYhat,'color',colores(3,:));
                
                edges = [0, tt(isnan(vR))];
                for i=1:length(edges)-1
                    xx = (edges(i)+edges(i+1))/2;
                    yy = -0.05;
                    ht(i) = text(xx,yy,['trial ',num2str(i)],'horizontalalignment','center');
                end
                
                ylabel('S_{when} [a.u.]');
                xlabel('Time [s]');
                set(gca,'tickdir','out','ticklength',[0.01,0.01]);
                ylim([-0.1,1.1]);
                set(gca,'xtick',[0,0.25],'xticklabel',{'250ms',''});
                %set(ht,'FontSize',10);
                
                here = pwd;
                [filepath , name] = fileparts(here);
                fname = fullfile(filepath, 'Final Figures', 'Figure_S6');
                saveas(gca, [fname '.fig'])
                saveas(gca, [fname '.pdf'])
                saveas(gca, [fname '.png'])
                
            end
            
        end
        
    end
    
    %%
    disp('Saving...')
    global mLoc
    savefilename = fullfile(mLoc, 'when_decoder', 'data_by_neuronal_group', dataset);
    save(savefilename,'W','names');
    
    %% Save only weight vector for When decoder
    w = W(:, 1);
    save(fullfile(saveLoc, ['WhenD_w_S' num2str(sess)]), 'w')
    
    allWhenDs{sess} = w;
end


%% Save
save(fullfile(saveLoc, 'AUC_by_neuronal_group'), 'AUC')

%% Plot goodness of fit

figure; hold on
n = length(names);
hb = bar(nanmean(AUC),'base',0.5);
set(hb,'FaceColor',0.7*[1,1,1],'EdgeColor','k');
hold all
terrorbar(1 : n,nanmean(AUC),stderror(AUC),'color','k','marker','o','markerfacecolor','k','linestyle','none');
set(gca,'xtick',1:n,'xticklabel',names);
ylabel('AUC');




