
medA = load(fullfile(saveLoc, 'for_plots_mediation_norm_to_se'));

mediatedDimsA = [52 19 27]; dimTin = 27;% ramp, PC1, Tin ('for_plots_mediation_norm_to_se')
tt = medA.leverage(medD).per_session.self_mediated.tt;
tpt = find(round(tt, 3) == 0.40);

allMed = [];
allP = [];

for medD = mediatedDimsA
    
    disp(medA.leverage(medD).signal)
    
    for sess = 1 : 8
        rCho(sess) = medA.leverage(medD).per_session.self_mediated(sess).rho_choice(tpt);
        rCho_part(sess) = medA.leverage(medD).per_session.self_mediated(sess).rho_choice_partial(tpt);
        rCho_part_Tin(sess) = medA.leverage(medD).per_session.other_mediated(sess).rho_choice_partial(tpt);
        
        rRT(sess) = medA.leverage(medD).per_session.self_mediated(sess).rho_RT(tpt);
        rRT_part(sess) = medA.leverage(medD).per_session.self_mediated(sess).rho_RT_partial(tpt);
        rRT_part_Tin(sess) = medA.leverage(medD).per_session.other_mediated(sess).rho_RT_partial(tpt);
    end
    
    medRT = 1 - rRT_part.^2./rRT.^2;
    medRT(medRT<0) = 0;
    disp('Mediation of RT effect')
    disp([num2str(nanmean(medRT)) ' +- ' num2str(nanstd(medRT)/sqrt(8))])
    allMed.RT{medD} = medRT;
    [h,p,ci,stats] = ttest(medRT);
    disp(['p=' num2str(p)])
    allP.RT{medD} = p;
    medRT = 1 - rRT_part_Tin.^2./rRT.^2;
    medRT(medRT<0) = 0;
    disp('Mediation of RT effect by Tin')
    disp([num2str(nanmean(medRT)) ' +- ' num2str(nanstd(medRT)/sqrt(8))])
    allMed.RT_tin{medD} = medRT;
    [h,p,ci,stats] = ttest(medRT);
    disp(['p=' num2str(p)])
    allP.RT_tin{medD} = p;
    
    
    
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
    allMed.Cho{medD} = medCho;
    [h,p,ci,stats] = ttest(medCho);
    allP.Cho{medD} = p;
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
    allMed.Cho_tin{medD} = medCho;
    [h,p,ci,stats] = ttest(medCho);
    allP.Cho_tin{medD} = p;
    disp(['p=' num2str(p)])
end


% At 400
% ramp
% Mediation of RT effect
% 0.8455 +- 0.039675
% p=1.2616e-07
% Mediation of RT effect by Tin
% 0.51082 +- 0.09925
% p=0.0013287
% Mediation of choice effect
% 0.71912 +- 0.097773
% p=0.00015526
% Mediation of choice effect by Tin
% 0.46741 +- 0.15023
% p=0.017052
% PC1
% Mediation of RT effect
% 0.74803 +- 0.067791
% p=1.1142e-05
% Mediation of RT effect by Tin
% 0.57761 +- 0.099065
% p=0.00064315
% Mediation of choice effect
% 0.60402 +- 0.12416
% p=0.0018249
% Mediation of choice effect by Tin
% 0.39008 +- 0.11265
% p=0.010509
% TinC
% Mediation of RT effect
% 0.57258 +- 0.092687
% p=0.00045511
% Mediation of RT effect by Tin
% 0.60945 +- 0.0598
% p=1.8878e-05
% Mediation of choice effect
% 0.67868 +- 0.10918
% p=0.00043826
% Mediation of choice effect by Tin
% 0.65947 +- 0.097596
% p=0.00026329




