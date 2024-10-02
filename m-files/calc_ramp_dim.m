% calc_ramp_dim
% Adaptend from AZ 2022
% Compute weights for ramp direction
% Includes weights computed with and without lasso regularization and with
% and without Tin neurons included in the population

global par_sess

save_fig = 1;
save_weights_flag = 1;
plot_fig = 0;
for use_lasso_flag = 0 : 1
    if use_lasso_flag
        fig_filename = 'output_lasso.pdf';
    else
        fig_filename = 'output.pdf';
    end
    
    %%
    deltat_ms = 25; % time bin
    
    v_dataset = 1 : 8;
    
    for sess = 1 : length(v_dataset)
        
        disp(['Session ' num2str(sess) ': Computing ramp direction...'])
        
        close all
        
        %% Load data
        
        par = par_sess{sess};
        dataset = sess;
        
        d = get_data(deltat_ms,[],sess, par,'dt_rel_RT',0);
        
        %%
        [nneurons, ntimes, ntrials] = size(d.H);
        
        %% construct the independent variable for the regression
        dt = d.t(2)-d.t(1);
        tind = d.t>=0.199 & d.t<=nanmax(d.RT+dt-0.05);
        m = min(d.t(tind));
        N = ceil(max((d.RT-m-0.05)/dt));
        Y = nan(N, ntrials);
        
        regress_to_flag = 1;
        
        for i=1:ntrials
            switch regress_to_flag
                case 1 % integration
                    rt = d.RT(i);
                    n = ceil((rt-m-0.05)/dt);
                    Y(1:n+1,i) = [0:n]/n;
                case 2 % momentary evidence
                    rt = d.RT(i);
                    n = ceil((rt-m-0.05)/dt);
                    Y(1:n+1,i) = 1/rt;
                case 3 % random stepping
                    rt = d.RT(i);
                    n = ceil((rt-m-0.05)/dt);
                    r = randsample(1:n,1);
                    Y(1:n+1,i) = 0;
                    Y(r:n+1,i) = 1;
            end
            
        end
        
        E = d.H(:,tind,:);
        
        %% regress
        
        % sanity check
        if size(E,2)~=size(Y,1)
            error('wrong sizes');
        end
        
        tr_ind = d.choice==0; % Tin choices
        
        Y = Y(:,tr_ind);
        Y = Y(:);
        E = E(:,:,tr_ind);
        E = E(:,:);
        
        if use_lasso_flag
            Y = 2*(Y - 0.5); % zero-centered
        end
        
        
        %%
        I = ~all(isnan(E))' & ~isnan(Y);
        depvar = Y(I);
        indepvar = {'neurons', E(:,I)', 'ones', ones(sum(I),1)};
        
        if use_lasso_flag
            lambda = 0.005;
            x = cat(2,indepvar{2:2:end});
            
            % normalize
            %         xn = nanzscore(x);% normalize
            media = nanmean(x);
            stdev = nanstd(x);
            xn = (x - media)./stdev;
            xn(:,end) = x(:,end); % keep the b_0 unnormalized
            J = media==0; % no spikes
            xn(:,J) = 0;
            
            [beta,stats] = lassoglm(xn,depvar,'normal','lambda',lambda);
            %         beta = lasso(x,depvar,'lambda',lambda);
        else
            [beta,idx,stats,x,LRT] = f_regression(depvar,[],indepvar);
        end
        
        
        w = beta(1:end-1);
        b0 = beta(end);
        if save_weights_flag
            if use_lasso_flag
                %savefilename = fullfile(filepath, 'All sessions', 'ramp_weights_lasso', d.dataset);
                savePath = fullfile(saveLoc, 'ramp_weights_lasso');
                if exist(savePath) == 0; mkdir(savePath); end
                
                save(fullfile(savePath, d.dataset),'w','b0','use_lasso_flag','media','stdev');
                % save(fullfile(filepath, 'All sessions', ['Ramp_w_S' num2str(dataset)]),'w','b0','use_lasso_flag','media','stdev');
                save(fullfile(saveLoc, ['Ramp_w_S' num2str(sess)]),'w','b0','use_lasso_flag','media','stdev');
                
            else
                % savefilename = fullfile(filepath, 'All sessions', 'ramp_weights', d.dataset);
                savePath = fullfile(saveLoc, 'ramp_weights');
                if exist(savePath) == 0; mkdir(savePath); end
                
                save(fullfile(savePath, d.dataset),'w','b0','use_lasso_flag');
                % save(fullfile(filepath, 'All sessions', ['RampNoLasso_w_S' num2str(dataset)]),'w','b0','use_lasso_flag');
                save(fullfile(saveLoc, ['RampNoLasso_w_S' num2str(sess)]),'w','b0','use_lasso_flag');
            end
        end
        
        %% Plot the fit for Session 1 with lasso regression (Supp Fig 2)
        if sess == 1 && use_lasso_flag == 1
            
            if use_lasso_flag
                yhat = glmval(beta,xn,'identity','constant','off');
            else
                yhat = glmval(beta,x,'identity','constant','off');
            end
            
            ff = find(depvar==min(depvar));
            ntr_for_plot = 6;
            first_trial = 300;
            figure('Position', [432  572  707  313]); hold on
            colores = movshon_colors(3);
            for k=1:ntr_for_plot
                ind = ff(first_trial + k - 1):[ff(first_trial + k)-1];
                ind0 = ff(first_trial);
                plot(ind-ind0, depvar(ind),'k');
                hold all
                plot(ind-ind0, yhat(ind),'color',colores(3,:));
            end
            ylabel('Activation');
            hl = legend('Expectation from DDM','Reconstruction');
            set(hl,'location','best');
            xlabel('Time step');
            set(gca,'tickdir','out','ticklength',[0.01,0.01]);
            xlim([0, max(ind-ind0)]);
            here = pwd;
            [filepath , name] = fileparts(here);
            fname = fullfile(filepath, 'Final Figures', 'Figure_S2');
            saveas(gca, [fname '.fig'])
            saveas(gca, [fname '.pdf'])
            saveas(gca, [fname '.png'])
        end
        %%
        h = w'*d.H(:,:) + b0;
        h = reshape(h,[ntimes, ntrials]);
          
        %% knock-out analysis
        %% no Tin
        w = beta(1:end-1);
        w(d.unitIdxLIP_Tin) = 0;
        h_noTin = w'*d.H(:,:) + b0;
        h_noTin = reshape(h_noTin,[ntimes, ntrials]);
        
        %% only Tins
        w = beta(1:end-1);
        J = ~ismember(1:length(w),d.unitIdxLIP_Tin);
        w(J) = 0;
        h_onlyTin = w'*d.H(:,:) + b0;
        h_onlyTin = reshape(h_onlyTin,[ntimes, ntrials]);
        
        
        %% refitting without Tin neurons
        
        ind = 1:nneurons;
        ind(ismember(ind, d.unitIdxLIP_Tin)) = [];
        
        indepvar = {'neurons', E(ind,I)', 'ones', ones(sum(I),1)};
        
        if use_lasso_flag
            x = cat(2,indepvar{2:2:end});
            xn = nanzscore(x);
            [beta,stats] = lassoglm(x,depvar,'normal','lambda',lambda);
        else
            [beta,idx,stats,x,LRT] = f_regression(depvar,[],indepvar);
        end
                
        w = zeros(nneurons,1);
        w(ind) = beta(1:end-1);
        b0 = beta(end);
        h_refit_noTin = w'*d.H(:,:) + b0;
        h_refit_noTin = reshape(h_refit_noTin,[ntimes, ntrials]);
       
        % Save w without Tin neurons
        if save_weights_flag
            if use_lasso_flag
                savePath = fullfile(saveLoc, 'ramp_weights_lasso');
                
                save(fullfile(savePath, [d.dataset '_noTin']),'w','b0','use_lasso_flag','media','stdev');
                save(fullfile(saveLoc, ['Ramp_noTin_w_S' num2str(sess)]),'w','b0','use_lasso_flag','media','stdev');
                
            else
                savePath = fullfile(saveLoc, 'ramp_weights');
                
                save(fullfile(savePath, [d.dataset '_noTin']),'w','b0','use_lasso_flag');
                save(fullfile(saveLoc, ['RampNoLasso_noTin_w_S' num2str(sess)]),'w','b0','use_lasso_flag');
            end
        end
        
    end
end

