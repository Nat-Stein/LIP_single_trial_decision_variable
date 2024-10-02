% mediation_zeta_control

load(fullfile(saveLoc, 'sig_allSessions'),...
    'sigSL')
S = sigSL;
% Analysis parameters
max_coh = 0.1;
minRT = 0.67;
maxRT = 2;
start_t = 0;
end_t = 0.55;

session = S.session;
times = S.t;
dt = times(2) - times(1);
choice = S.choice;
RT = S.rt;
coh = S.sig_coh;

str =  {'PC1',...
    'TinC',...
    'ramp'};

% Run mediation only on one time point!


plotFigs = 0;

% smooth in window
for i = 1 : length(str)
    sm = round(0.05/dt);
    h = ones(sm,1)/sm;
    S.(str{i}) = conv2(1, h, S.(str{i}), 'same');
end

%% time subset
tind = findclose(times, -0.1:0.01:0.75);
for i=1:length(str)
    aux = S.(str{i});
    S.(str{i}) = aux(:,tind);
end
times = times(tind);


%% Run mediation itself

warning('off','all')

allCoh = unique(coh)';
flags.norm_to_se = 1;
nsessions =  8;
clear lev_singleC
for cc = 1 : length(allCoh)
    sig_coh = allCoh(cc);
    disp(['Coherence ' num2str(sig_coh)])
    clear leverage
    count = 0;
    for iMethod=1:length(str)
        sProj = S.(str{iMethod});
        disp(['Method ' num2str(iMethod) '/' num2str(length(str))]);
        for iMediator=1:length(str)
            
            if iMediator~=iMethod
                
                tind = findclose(times, 0.55);
                sProj_mediator = sProj(:,tind);
                sOther_mediator = S.(str{iMediator})(:,tind);
                
                count = count + 1;
                
                %% per session and average
                do_plot = 0;
                clear out out_other
                i = 0; incl_sess = [];
                for j=1:nsessions
                    I = session==j & coh == sig_coh & RT > minRT; % only low coh
                    numTrl(cc, j) = sum(I);
                    if ~[all(isnan(to_vec(sProj(I,:)))) || all(isnan(sProj_mediator(I))) || all(isnan(sOther_mediator(I)))] & sum(I) > 5 & abs(sig_coh) < 0.2
                        i = i+1;
                        incl_sess = [incl_sess j];
                        [~,out(i)] = corr_with_RT_choice(choice(I)==0, RT(I), coh(I), times, ...
                            sProj(I,:), sProj_mediator(I), minRT, maxRT, start_t, end_t, do_plot,flags);
                        [~,out_other(i)] = corr_with_RT_choice(choice(I)==0, RT(I), coh(I), times, ...
                            sProj(I,:), sOther_mediator(I), minRT,maxRT, start_t, end_t,do_plot,flags);
                    end
                end
                
                %% save
                if i > 0
                    leverage(count).signal = str{iMethod};
                    leverage(count).mediator = str{iMediator};
                    leverage(count).incl_sess = incl_sess;
                    
                    leverage(count).t = out(1).tt;
                    
                    
                    x = cat(2, out.rho_RT); 
                    [leverage(count).RT.unmediated.m , leverage(count).RT.unmediated.se] = averageCorrelation(x);
                    
                    x = cat(2, out.rho_RT_partial); 
                    [leverage(count).RT.mediated_self.m, leverage(count).RT.mediated_self.se] = averageCorrelation(x);
                    
                    x = cat(2, out_other.rho_RT_partial); 
                    [leverage(count).RT.mediated_other.m, leverage(count).RT.mediated_other.se] = averageCorrelation(x);
                    
                    
                    leverage(count).choice.unmediated.m = nanmean(cat(2,out.rho_choice),2);
                    leverage(count).choice.unmediated.se = stderror(cat(2,out.rho_choice),2);
                    leverage(count).choice.mediated_self.m = nanmean(cat(2,out.rho_choice_partial),2);
                    leverage(count).choice.mediated_self.se = stderror(cat(2,out.rho_choice_partial),2);
                    leverage(count).choice.mediated_other.m = nanmean(cat(2,out_other.rho_choice_partial),2);
                    leverage(count).choice.mediated_other.se = stderror(cat(2,out_other.rho_choice_partial),2);
                    
                    leverage(count).per_session.self_mediated = out;
                    leverage(count).per_session.other_mediated = out_other;
                end
                
                
                
            end
        end
    end
    if abs(sig_coh) < 0.2
        save(fullfile(saveLoc, ['leverage_cc' num2str(cc)]), 'leverage', 'start_t', 'end_t')

        lev_singleC{cc}.leverage = leverage;
        save(fullfile(saveLoc, 'lev_singleC'), 'lev_singleC', 'start_t', 'end_t', 'numTrl')
    end
end

%% Compute stats on leverage
K = lev_singleC;
tt = K{cc}.leverage(1).per_session.self_mediated(1).tt;
tpt = find(round(tt, 3) == 0.40);

count = 1;

clear pvals
allCoh = unique(coh)';
numtrl = [];
for count = 1 : 6
    
    i = 0;
    for cc = [3 : 9]
        clear p x y
        sig_coh = allCoh(cc);
        i = i + 1;
        for sess = 1 : 8
            x(sess) = K{cc}.leverage(count).per_session.self_mediated(sess).rho_RT(tpt);
            y(sess) = K{cc}.leverage(count).per_session.self_mediated(sess).rho_choice(tpt);
            
            for ch = 0 : 1
                I = session==sess & coh == sig_coh & RT > minRT & choice == ch;
                numtrl(i, sess, ch+1) = sum(I);
            end
            
        end
        x = x(~isnan(x));
        [avg_corr, sem_corr, z_scores, p] = averageCorrelation(x);
        pvals.RT(count, i) = p;
        
        y = y(~isnan(y));
        [h,p,ci,stats] = ttest(y);
        pvals.choice(count, i) = p;
        
    end
end

%% ---------------------------------------------------------------------------
% Compute mediation stats

K = lev_singleC;

mediatedDimsA = [6 1 3]; dimTin = 3; % ramp, PC1, Tin

tt = K{cc}.leverage(1).per_session.self_mediated(1).tt;
tpt = find(round(tt, 3) == 0.40);

allMed = [];
allP = [];
for cc = 3 : 9
    for medD = mediatedDimsA
        
        disp([K{cc}.leverage(medD).signal ', coh=' num2str(allCoh(cc))])
        
        for sess = 1 : 8
            rCho(sess) = K{cc}.leverage(medD).per_session.self_mediated(sess).rho_choice(tpt);
            rCho_part(sess) = K{cc}.leverage(medD).per_session.self_mediated(sess).rho_choice_partial(tpt);
            rCho_part_Tin(sess) = K{cc}.leverage(medD).per_session.other_mediated(sess).rho_choice_partial(tpt);
            
            rRT(sess) = K{cc}.leverage(medD).per_session.self_mediated(sess).rho_RT(tpt);
            rRT_part(sess) = K{cc}.leverage(medD).per_session.self_mediated(sess).rho_RT_partial(tpt);
            rRT_part_Tin(sess) = K{cc}.leverage(medD).per_session.other_mediated(sess).rho_RT_partial(tpt);
        end
        
        medRT = 1 - rRT_part.^2./rRT.^2;
        medRT(medRT<0) = 0;
        disp('Mediation of RT effect')
        disp([num2str(nanmean(medRT)) ' +- ' num2str(nanstd(medRT)/sqrt(8))])
        allMed.RT{medD, cc} = medRT;
        [h,p,ci,stats] = ttest(medRT);
        disp(['p=' num2str(p)])
        allP.RT(medD, cc) = p;
        medRT = 1 - rRT_part_Tin.^2./rRT.^2;
        medRT(medRT<0) = 0;
        disp('Mediation of RT effect by Tin')
        disp([num2str(nanmean(medRT)) ' +- ' num2str(nanstd(medRT)/sqrt(8))])
        allMed.RT_tin{medD, cc} = medRT;
        [h,p,ci,stats] = ttest(medRT);
        disp(['p=' num2str(p)])
        allP.RT_tin(medD, cc) = p;
        
        
        for i = 1 : length(rCho)
            if rCho_part(i) <= 0 & rCho(i) > 0
                medCho(i) = 1;
            elseif rCho(i) <= 0
                medCho(i) = nan;
            else
                medCho(i) = 1 - rCho_part(i)/rCho(i);
            end
        end
        medCho(medCho<0) = 0;
        disp('Mediation of choice effect')
        disp([num2str(nanmean(medCho)) ' +- ' num2str(nanstd(medCho)/sqrt(8))])
        allMed.Cho{medD, cc} = medCho;
        [h,p,ci,stats] = ttest(medCho);
        allP.Cho(medD, cc) = p;
        disp(['p=' num2str(p)])
        
        for i = 1 : length(rCho)
            if rCho_part_Tin(i) <= 0 & rCho(i) > 0
                medCho(i) = 1;
            elseif rCho(i) <= 0
                medCho(i) = nan;
            else
                medCho(i) = 1 - rCho_part_Tin(i)/rCho(i);
            end
        end
        medCho(medCho<0) = 0;
        disp('Mediation of choice effect by Tin')
        disp([num2str(nanmean(medCho)) ' +- ' num2str(nanstd(medCho)/sqrt(8))])
        allMed.Cho_tin{medD, cc} = medCho;
        [h,p,ci,stats] = ttest(medCho);
        allP.Cho_tin(medD, cc) = p;
        disp(['p=' num2str(p)])
    end
end


max(allP.RT)
max(allP.RT_tin)

max(allP.Cho)
max(allP.Cho_tin)


% Just do: are r-values and betas different between leverage and mediation






