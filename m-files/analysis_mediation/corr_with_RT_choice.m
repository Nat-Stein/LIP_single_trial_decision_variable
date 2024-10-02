function [p,out] = corr_with_RT_choice(choice, RT, coh, times, dv, dv_mediator, minRT,maxRT, start_t, end_t,do_plot,flags)


if do_plot ==0
    p = [];
end

% 
norm_to_se_flag = flags.norm_to_se; % testing

%
filt = RT>minRT & RT<maxRT;
ucoh = nanunique(coh);

% residuals
dv_res = nan(size(dv));
dv_res_choice = nan(size(dv));
dv_res_med_choice = nan(size(dv_mediator));
dv_res_med = nan(size(dv_mediator));
uchoice = [0,1];
RT_res = nan(size(RT));
for i=1:length(ucoh)
    I = coh==ucoh(i) & filt;
    if sum(I)>0
        dv_res_choice(I,:) = dv(I,:) - nanmean(dv(I,:));
        dv_res_med_choice(I) = dv_mediator(I) - nanmean(dv_mediator(I));
    end
    for j=1:length(uchoice)
        I = choice==uchoice(j) & coh==ucoh(i) & filt;
        if sum(I)>0
            dv_res(I,:) = dv(I,:) - nanmean(dv(I,:));
            RT_res(I) = RT(I) - nanmean(RT(I));
            dv_res_med(I) = dv_mediator(I) - nanmean(dv_mediator(I));
        end
    end
end

%%

tind = find(times>=start_t & times<=end_t);
tt = times(tind);
% ind_mediator = findclose(times, 0.55);
out.tt = tt;



%% parcorr with RT w/wo mediation

idx = filt & choice==1;
if sum(idx) > 2
    rho_RT = corr(dv_res(idx,tind),RT_res(idx),'rows','complete');
    rho_RT_partial = partialcorr(dv_res(idx,tind),RT_res(idx),...
        dv_res_med(idx),'rows','complete');
    
    if do_plot
        p = publish_plot(2,1);
        p.next();
        
        
        plot(tt, rho_RT)
        hold all
        plot(tt, rho_RT_partial)
        ylabel('\rho with RT');
        legend('unmediated','mediated');
    end
    
    out.rho_RT = rho_RT;
    out.rho_RT_partial = rho_RT_partial;
else
    out.rho_RT = nan(length(tind), 1);
    out.rho_RT_partial = nan(length(tind), 1);
end
%% beta with choice with/without mediation
idx = filt;


rho_choice = nan(length(tind),1);
rho_choice_partial = nan(length(tind),1);
norm_choice = nan(length(tind),1);
norm_choice_partial = nan(length(tind),1);
for i=1:length(tind)
    indepvar = {'dv_res_choice',dv_res_choice(idx,tind(i)), ...
        'ones',ones(sum(idx),1)};
    [beta,indices,stats,x,LRT] = ...
        f_regression(choice(idx),[],indepvar);
    
    
    if norm_to_se_flag
        norm = stats.se(indices.dv_res_choice);
        rho_choice(i) = beta(indices.dv_res_choice)./norm;
    else
        rho_choice(i) = beta(indices.dv_res_choice);
    end
    norm_choice(i) = norm;
%     indepvar = {'dv_res_choice',dv_res_choice(idx,tind(i)), ...
%         'mediator',dv_res_med_choice(idx) ,'ones',ones(sum(idx),1)}; % to do: try this with one bias term per coh
    indepvar = {'dv_res_choice',dv_res_choice(idx,tind(i)), ...
        'mediator',dv_res_med_choice(idx) ,'ones',adummyvar(coh(idx))}; 

    [beta,indices,stats,x,LRT] = ...
        f_regression(choice(idx),[],indepvar);
    
    if norm_to_se_flag
        rho_choice_partial(i) = beta(indices.dv_res_choice)./norm; % use the norm of the unmediated
    else
        rho_choice_partial(i) = beta(indices.dv_res_choice);
    end
    norm_choice_partial(i) = norm;
end

out.rho_choice = rho_choice;
out.rho_choice_partial = rho_choice_partial;
out.norm_choice = norm_choice;
out.norm_choice_partial = norm_choice_partial;

if do_plot
    p.next();
    plot(tt,rho_choice, tt, rho_choice_partial);



    ylabel('\beta with choice');
    xlabel('Time [s]');
end



% plot(tt(1:5:end), beta(idx.dv_res_choice));

%%
