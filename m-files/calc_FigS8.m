medA = load(fullfile(saveLoc, 'for_plots_mediation_norm_to_se'));

uni_s = unique({medA.leverage.signal});
labels = replace(uni_s,'w_','');

m = {medA.leverage.mediator};
s = {medA.leverage.signal};

X = nan(length(uni_s));
medCho_other_ForRatio = nan(2, 8, 8); medCho_raw_ForRatio = medCho_other_ForRatio; medRT_raw_ForRatio = medCho_other_ForRatio;
clear X_ratioS Xrt_ratioS X_ratioO Xrt_ratioO medCho_other_ForRatio medRT_other_ForRatio
for i = 1 : length(uni_s)
    for j = 1 : length(uni_s)
        
        
        if i == j
            I = find(ismember(s,uni_s(i)));
            L = medA.leverage(I(1));
            t = single(L.t);
            tind1 = findclose(t,0.4);
            
            numSess = length(L.per_session.self_mediated);
            % Choice
            clear raw m_self m_other medCho_self_stats
            for sess = 1 : numSess
                % Choice
                if length(L.per_session.self_mediated(sess).rho_choice) > 0
                    raw(sess) = L.per_session.self_mediated(sess).rho_choice(tind1);
                    m_self(sess) = L.per_session.self_mediated(sess).rho_choice_partial(tind1);
                else
                    raw(sess) = nan;
                    m_self(sess) = nan;
                end
            end
            
            % Self
            forRatio = nan(1, 8);
            val_tpt = raw > 0 & m_self > 0 & raw > m_self;
            medCho_self_stats(val_tpt) = 1 - (m_self(val_tpt) ./ raw(val_tpt));
            forRatio(val_tpt) = m_self(val_tpt);
            inval_tpt = raw > 0 & m_self > 0 & raw < m_self;
            medCho_self_stats(inval_tpt) = 0;
            forRatio(inval_tpt) = raw(inval_tpt);
            inval_tpt = raw > 0 & m_self < 0;
            medCho_self_stats(inval_tpt) = 1;
            forRatio(inval_tpt) = 0;
            inval_tpt = raw < 0;
            medCho_self_stats(inval_tpt) = nan;
            forRatio(inval_tpt) = nan;
            
            mediation = nanmean(medCho_self_stats);
            medCho_other_ForRatio(i, i, :) = forRatio;
            medCho_raw_ForRatio(i, i, 1 : numSess) = raw;
            
            % RT
            clear raw m_self m_other medRT_self medRT_self_ForRatio
            for sess = 1 : numSess
                % RT
                if length(L.per_session.self_mediated(sess).rho_RT) > 0
                    raw(sess) = L.per_session.self_mediated(sess).rho_RT(tind1);
                    m_self(sess) = L.per_session.self_mediated(sess).rho_RT_partial(tind1);
                else
                    raw(sess) = nan;
                    m_self(sess) = nan; 
                end
            end
            forRatio = nan(1, 8);
            val_tpt = raw < 0 & m_self <0 & raw < m_self;
            medRT_self(val_tpt) = 1 - ((m_self(val_tpt)).^2 ./ (raw(val_tpt)).^2);
            forRatio(val_tpt) = m_self(val_tpt);
            inval_tpt = raw < 0 & m_self <0 & raw > m_self;
            medRT_self(inval_tpt) = 0;
            forRatio(inval_tpt) = raw(inval_tpt);
            inval_tpt = raw < 0 & m_self >0;
            medRT_self(inval_tpt) = 1;
            forRatio(inval_tpt) = 0;
            inval_tpt = raw > 0;
            medRT_self(inval_tpt) = nan;
            forRatio(inval_tpt) = nan;
            
            mediationRT = nanmean(medRT_self);
            medRT_other_ForRatio(i, i, :) = forRatio;
            medRT_raw_ForRatio(i, i, 1 : numSess) = raw;
%             X_ratio(i,j) = 1;
%             Xrt_ratio(i,j) = 1;
        else
            I = ismember(s,uni_s(i)) & ismember(m,uni_s(j));
            L = medA.leverage(I);
            t = single(L.t);
            tind1 = findclose(t,0.4);
            
            numSess = length(L.per_session.other_mediated);
            
            clear raw m_self m_other medCho_other_stats
            for sess = 1 : numSess
                raw1 = L.per_session.self_mediated(sess).rho_choice;
                m_other1 = L.per_session.other_mediated(sess).rho_choice_partial;
                m_self1 = L.per_session.self_mediated(sess).rho_choice_partial;
                if length(raw1) > 0 & length(m_other1) > 0 & length(m_self1) > 0
                    raw(sess) = raw1(tind1);
                    m_other(sess) = m_other1(tind1);
                    m_self(sess) = m_self1(tind1);
                else
                    raw(sess) = nan;
                    m_other(sess) = nan;
                    m_self(sess) = nan;
                end
                
            end
            
            % Other
            forRatio = nan(1, 8);
            val_tpt = raw > 0 & m_other > 0 & raw > m_other;
            medCho_other_stats(val_tpt) = 1 - (m_other(val_tpt) ./ raw(val_tpt));
            forRatio(val_tpt) = m_other(val_tpt);
            % medCho_other_ForRatio(val_tpt) = (m_other(val_tpt) - raw(val_tpt)) ./ (m_self(val_tpt) - raw(val_tpt));
            inval_tpt = raw > 0 & m_other > 0 & raw < m_other;
            medCho_other_stats(inval_tpt) = 0;
            forRatio(inval_tpt) = raw(inval_tpt);
            inval_tpt = raw > 0 & m_other < 0;
            medCho_other_stats(inval_tpt) = 1;
            forRatio(inval_tpt) = 0;
            inval_tpt = raw < 0;
            medCho_other_stats(inval_tpt) = nan;
            forRatio(inval_tpt) = nan;
            
            mediation = nanmean(medCho_other_stats);
            
            % Mediation ratio
            medCho_other_ForRatio(i, j, :) = forRatio; %(m_other - raw) ./ (m_self - raw);
            medCho_raw_ForRatio(i, j, 1 : numSess) = raw;
            % ---------------------------------------------------------------------------
            % RT
            clear raw m_self m_other medRT_other 
            for sess = 1 : length(L.per_session.other_mediated)
                raw1 = L.per_session.self_mediated(sess).rho_RT;
                m_other1 = L.per_session.other_mediated(sess).rho_RT_partial;
                m_self1 = L.per_session.self_mediated(sess).rho_RT_partial;
                if length(raw1) > 0 && length(m_other1) > 0 && length(m_self1) > 0
                    raw(sess) = raw1(tind1);
                    m_other(sess) = m_other1(tind1);
                    m_self(sess) = m_self1(tind1);
                else
                    raw(sess) = nan;
                    m_other(sess) = nan;
                    m_self(sess) = nan;
                end
            end
            % Other - should be different equation!
            forRatio = nan(1, 8);
            val_tpt = raw < 0 & m_other <0 & raw < m_other;
            medRT_other(val_tpt) = 1 - ((m_other(val_tpt)).^2 ./ (raw(val_tpt)).^2);
            forRatio(val_tpt) = m_other(val_tpt);
            inval_tpt = raw < 0 & m_other <0 & raw > m_other;
            medRT_other(inval_tpt) = 0;
            forRatio(inval_tpt) = raw(inval_tpt);
            inval_tpt = raw < 0 & m_other >0;
            medRT_other(inval_tpt) = 1;
            forRatio(inval_tpt) = 0;
            inval_tpt = raw > 0;
            medRT_other(inval_tpt) = nan;
            forRatio(inval_tpt) = nan;
            
            mediationRT = nanmean(medRT_other);
            
            % Mediation ratio
            medRT_other_ForRatio(i, j, :) = forRatio; %(m_other - raw) ./ (m_self - raw);
            medRT_raw_ForRatio(i, j, 1 : numSess) = raw;
        end
        X(i,j) = mediation;
        Xrt(i,j) = mediationRT;
        
%         X_ratio(i,j) = nanmean(medCho_other_ForRatio(i, j, :));
%         Xrt_ratio(i,j) = nanmean(medRT_other_ForRatio(i, j, :));
        
    end
end
% YYYYYYYYYYYYYYYYYY - something is still missing here
clear X_ratio Xrt_ratio
for i = 1 : size(X, 1)
    for j = 1 : size(X, 2)
        if i == j
            X_ratio(i,j) = 1;
            Xrt_ratio(i,j) = 1;
        else
            X_ratio(i,j) = X(i, j) / X(i, i);
            Xrt_ratio(i,j) = Xrt(i, j) / Xrt(i, i);
        end
    end
end

save(fullfile(saveLoc, 'xMediationMatrices'), 'X', 'Xrt', 'X_ratio', 'Xrt_ratio')
% save('E:\Matalb analyses\xMediationMatrices_230522', 'X', 'Xrt', 'X_ratio', 'Xrt_ratio')


