addpath('../generic/');
addpath(genpath('../az_matlab_files/'));

%%

ndatasets = 8;

names = {'all','Tin','Tout','Min_C','Min_i','allButTin','allButTinOut'};
AUC = nan(ndatasets, length(names));

for idataset = 1:ndatasets

    dt_ms = 25; % 100
    %     dt_ms = 50;
    dataset = idataset;
    dat = get_data(dt_ms,[],dataset,'dt_rel_RT',0.05);

    struct2vars(dat);

    idx_include = choice==0;
    %     idx_include = choice==1;

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

    %% regre

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
            AUC(idataset,k) = auc;

            %%
            W(idx_neurons{k},k) = B;
            %             Weights{k} = B; % regression weights

            %%
            if k==1 % all neurons

                ntr_for_plot = 8;
                first_trial = 350;
                tr = unique(trial_id);
                tr_ind = findclose(tr, first_trial);
                tr_ind = tr(tr_ind:tr_ind + ntr_for_plot - 1);

                p = publish_plot(1,1);
                set(gcf,'Position',[432  572  707  313]);
                colores = movshon_colors(3);
                
                
                vR = [];
                vYhat = [];
                for i=1:length(tr_ind)
                    ind = trial_id ==tr_ind(i) & T>0.2;
                    vR = [vR; R(ind); NaN];
                    vYhat = [vYhat; yhat(ind); NaN];
                end
                
                
%                 for i=1:length(edges)-1
%                     ha = area([edges(i),edges(i+1)],[1,1]);
%                     set(ha,'edgecolor','w','facecolor',0.8 *[1,1,1] + iseven(i)*0.1*[1,1,1] );
%                 end
                

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
                %             set(hl,'location','best');
                xlabel('Time [s]');
                p.format('LineWidthPlot',2,'FontSize',14);
                set(gca,'tickdir','out','ticklength',[0.01,0.01]);

                

                %     save traces_for_plotting depvar yhat
                ylim([-0.1,1.1]);
                p.text_draw_fig(dat.dataset);

                p.offset_axes;
                p.visible_limits_axis(1,1,[0,0.25]);
                set(gca,'xtick',[0,0.25],'xticklabel',{'250ms',''});
                set(ht,'FontSize',10);
                p.append_to_pdf('fig_train',k==1 & idataset==1,1);

                
%                 close all

            end

            %% plot the dep var and the fit
            if 1
                if k==1 % all neurons
    
    
                    ntr_for_plot = 14;
    
                    first_trial = 350;
    
                    p = publish_plot(ntr_for_plot,2);
                    set(gcf,'Position',[584  275  682  726]);
                    for l=1:2
    
                        tr = unique(trial_id);
                        if l==1
                            tr(~iseven(tr)) = [];
                        else
                            tr(iseven(tr)) = [];
                        end
                        tr_ind = findclose(tr,first_trial);
                        tr_ind = tr(tr_ind:tr_ind + ntr_for_plot - 1);
    
                        colores = movshon_colors(3);
    
                        for itr=1:ntr_for_plot
    
                            p.next();
                            vR = [];
                            vYhat = [];
    
                            ind = trial_id ==tr_ind(itr) & T>=-0.2;
    
                            tt = dt * [1:sum(ind)] - dt;
                            stairs(T(ind), R(ind),'k');
                            hold all
                            plot(T(ind), yhat(ind),'color',colores(3,:));
                            
                            hold on
                            plot([RT(tr_ind(itr)), RT(tr_ind(itr))], [-0.1,1.1],'r');
                        end
    
                        %                 ylabel('Activation');
                        %             set(hl,'location','best');
                        %                 xlabel('Time [s]');
                        
    
                    end
    
                    set(p.h_ax,'xlim',[0.2,2]);
    
    %                 same_xscale(p.h_ax);
    
                    p.format('LineWidthPlot',1,'FontSize',9);
                    set(gca,'tickdir','out','ticklength',[0.01,0.01]);
    
                    %     save traces_for_plotting depvar yhat
                    set(p.h_ax,'ylim',[-0.1,1.1]);
                    p.text_draw_fig(1,'Train');
                    p.text_draw_fig(2,'Test');
                    p.text_draw_fig(dat.dataset);
                    p.unlabel_center_plots();
    
    
                    p.append_to_pdf('fig_train',0,1);
                    close all
                end
            end
        end

    end

    %%
    savefilename = fullfile('./data_by_neuronal_group',dataset);
    save(savefilename,'W','names');

end


%%
save AUC_by_neuronal_group AUC

%%
p = publish_plot(1,1);
n = length(names);
hb = bar(nanmean(AUC),'base',0.5);
set(hb,'FaceColor',0.7*[1,1,1],'EdgeColor','k');
hold all
terrorbar(1:n,nanmean(AUC),stderror(AUC),'color','k','marker','o','markerfacecolor','k','linestyle','none');
set(gca,'xtick',1:n,'xticklabel',names);
ylabel('AUC');
p.format();
p.append_to_pdf('fig_regre_to_when_by_neuronal_group',1,1);



