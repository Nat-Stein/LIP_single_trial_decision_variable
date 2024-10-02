% run_mediation_all2all
% by Ariel Zylberberg

load(fullfile(saveLoc, 'sig_allSessions'),...
    'sigSL')
S = sigSL;
% Analysis parameters
max_coh = 0.1;
minRT = 0.67;
maxRT = 2;
start_t = 0;
end_t = 0.6;

session = S.session;
times = S.t;
dt = times(2) - times(1);
choice = S.choice;
RT = S.rt;
coh = S.sig_coh;

S = sigSL;
S.MinC_minus_I = S.MinC - S.MinI;
S.Int_MinC_minus_I = cumsum(S.MinC - S.MinI, 2);

BL = find(times >= -0.1 & times < 0);
BL_minCI = nanmean(S.MinC_minus_I(:, BL), 2);
dat = S.MinC_minus_I - repmat(BL_minCI, 1, size(S.MinC_minus_I, 2));
dat(:, find(par.tsl<0)) = 0;

intMinCI_BL_orig = cumsum(dat, 2);
% intMinCI_BL_orig = S.Int_MinC_minus_I_BL;
A = nan(size(intMinCI_BL_orig));
doNorm = 1;
fw = ones(1, 50)/50;
tsl = find(par.tsl >= 0 & par.tsl <=0.6);
% Normalize amplitude of Int_MinC_minus_I_BL
if doNorm == 1
    for sess = 1 : 8
        
        trials = find(sigSL.session == sess);
        thisS = intMinCI_BL_orig(trials, :);
        cho = sigSL.choice(trials);
        
        mm = [];
        for ch = 0 : 1
            c0 = find(cho == ch);
            add = conv(squeeze(nanmean(thisS(c0, :), 1)), fw, 'same');
            mm = [mm add(1, tsl)];
            %add = conv(squeeze(nanmean(newRL(c0, :), 1)), fw, 'same');
            %mm = [mm add(1, trl)];
        end
        maxM = max(mm);
        minM = min(mm);
        exc = maxM - minM;
        
        newSL = (thisS - minM) ./ exc;
        
        A(trials, :) = newSL;
    end
end
S.Int_MinC_minus_I_BL = intMinCI_BL_orig;
S.Int_MinC_minus_I_BL_norm = A;

S.minReg = sigSL.minReg;
BL = find(times >= -0.1 & times < 0);
BL_minReg = nanmean(S.minReg(:, BL), 2);
dat = S.minReg - repmat(BL_minReg, 1, size(S.minReg, 2));
dat(:, find(par.tsl<0)) = 0;
S.Int_minReg = cumsum(dat, 2);

str =  {'MinC',...
    'MinI',...
    'PC1',...
    'TinC',...
    'whatD',...
    'whenC',...
    'ramp',...
    'MinC_minus_I',...
    'Int_MinC_minus_I',...
    'Int_MinC_minus_I_BL',...
    'Int_MinC_minus_I_BL_norm',...
    'Int_minReg'};

% S.PC1min = sigSL.PC1min;
% str =  {'MinC',...
%     'MinI',...
%     'PC1',...
%     'TinC',...
%     'whatD',...
%     'whenC',...
%     'ramp',...
%     'MinC_minus_I',...
%     'Int_MinC_minus_I',...
%     'Int_MinC_minus_I_BL',...
%     'PC1min'};


% 
% 
% str =  {'MinC',...
%     'MinI',...
%     'PC1',...
%     'TinC',...
%     'whatD',...
%     'whenC',...
%     'ramp',...
%     'MinC_minus_I',...
%     'Int_MinC_minus_I',...
%     'Int_MinC_minus_I_BL',...
%     'minReg',...
%     'Int_minReg'};

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


%%
warning('off','all')

flags.norm_to_se = 1;
nsessions =  8;

count = 0; % count = 111;
for iMethod=1:length(str)
    sProj = S.(str{iMethod});
    disp(['Method ' num2str(iMethod) '/' num2str(length(str))]);
    for iMediator=1:length(str)
        
        if iMediator~=iMethod
            
            tind = findclose(times, 0.55);
            sProj_mediator = sProj(:,tind);
            sOther_mediator = S.(str{iMediator})(:,tind);
            
            count = count + 1;
            
            %% all sessions together
            
            
            
            %% per session and average
            do_plot = 0;
            clear out out_other
            i = 0;
            for j=1:nsessions
                I = session==j & abs(coh)<max_coh; % only low coh
                if ~[all(isnan(to_vec(sProj(I,:)))) || all(isnan(sProj_mediator(I))) || all(isnan(sOther_mediator(I)))]
                    i = i+1;
                    [~,out(i)] = corr_with_RT_choice(choice(I)==0, RT(I), coh(I), times, ...
                        sProj(I,:), sProj_mediator(I), minRT,maxRT, start_t, end_t, do_plot,flags);
                    [~,out_other(i)] = corr_with_RT_choice(choice(I)==0, RT(I), coh(I), times, ...
                        sProj(I,:), sOther_mediator(I), minRT,maxRT, start_t, end_t,do_plot,flags);
                end
            end
            
            %% save
            
            leverage(count).signal = str{iMethod};
            leverage(count).mediator = str{iMediator};
            
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

% save(fullfile(saveLoc, 'for_plots_mediation_norm_to_se_MinCI_BL'), 'leverage')

save(fullfile(saveLoc, 'for_plots_mediation_norm_to_se_MinCI_BL'), 'leverage')














