
useDims = unique([mainDim compDims]); 
sc= 1;
for sess = 1 : 8
    for dim = useDims
        par = par_sess{sess};
        cohs = unique(par.dataMat(:, par.idx.sig_coh));
        % Create residuals
        for sc = 1% : size(allDim_sl{dim}, 1)
            dMat = allDim_sl{dim}{sc, sess};
            if size(dMat, 1) > 0
                rMat = nan(size(dMat));
                for c = 1 : length(cohs)
                    coh = cohs(c);
                    trls = find(par.dataMat(:, par.idx.sig_coh) == coh);
                    
                    rMat(trls, :) = dMat(trls, :) - repmat(nanmean(dMat(trls, :), 1), length(trls), 1);
                end
                rCell{dim, sess, sc} = rMat;
                rCellconv{dim, sess, sc} = conv2(1, filtW, rMat);
                rCellconv50{dim, sess, sc} = conv2(1, filtW50, rMat);
                rCellconv50b{dim, sess, sc} = conv2(1, filtW50b, rMat);
            end
        end
    end
end

rtAdd = 0.1; %used to be 0.08; now 0.1 again tp match other analyses % same as for RT and choice correlatinos, used to be 0.1 for otehr analyses. 
for sess = 1 : 8
    for dim1 = mainDim
        sc1 = 1; sn1 = 1; %if dim1 == 20 && sess == 2; sc1 = 2; sn1 = -1;  end
        for dim2 = compDims
            sc2 = 1; sn2 = 1; %if dim2 == 20 && sess == 2; sc2 = 2; sn2 = -1; end
            par = par_sess{sess};
            for sc2 = 1 %: size(allDim_sl{dim2}, 1)
                sn2 = 1; %if sess == 2 && dim == 20 && sc2 == 2; sn2 = -1; end
                if size(rCell{dim2, sess, sc2}, 1) > 0 &&  size(rCell{dim1, sess, sc1}, 1)
                    resCorr_wiTrial{dim1, dim2, sess, sc2} = nan(1, size(rCell{dim2, sess, sc2}, 1));
                    resCorr_wiTrialConv{dim1, dim2, sess, sc2} = nan(1, size(rCell{dim2, sess, sc2}, 1));
                    resCorr_wiTrialConv50bStep{dim1, dim2, sess, sc2} = nan(1, size(rCell{dim2, sess, sc2}, 1));
                    resCorr_wiTrialConv50{dim1, dim2, sess, sc2} = nan(1, size(rCell{dim2, sess, sc2}, 1));
                    resCorr_wiTrialConv50b{dim1, dim2, sess, sc2} = nan(1, size(rCell{dim2, sess, sc2}, 1));
                    for trl = 1 : size(rCell{dim2, sess, sc2}, 1)
                        tpts = find(~isnan(rCell{dim1, sess, sc1}(trl, :)) & ~isnan(rCell{dim2, sess, sc2}(trl, :)) &...
                            par.tsl >=0.2 & par.tsl <= par.dataMat(trl, par.idx.rt)/1000 - rtAdd); % Is nan as or 0.9s before saccade
                        
                        if length(tpts) > 5
                            [R, P] = corrcoef(sn1*rCell{dim1, sess, sc1}(trl, tpts), sn2*rCell{dim2, sess, sc2}(trl, tpts), 'rows', 'complete');
                            resCorr_wiTrial{dim1, dim2, sess, sc2}(trl) = R(1, 2);
                            
                            [R, P] = corrcoef(sn1*rCellconv{dim1, sess, sc1}(trl, tpts), sn2*rCellconv{dim2, sess, sc2}(trl, tpts), 'rows', 'complete');
                            resCorr_wiTrialConv{dim1, dim2, sess, sc2}(trl) = R(1, 2);
                            
                            [R, P] = corrcoef(sn1*rCellconv50{dim1, sess, sc1}(trl, tpts), sn2*rCellconv50{dim2, sess, sc2}(trl, tpts), 'rows', 'complete');
                            resCorr_wiTrialConv50{dim1, dim2, sess, sc2}(trl) = R(1, 2);
                            
                            [R, P] = corrcoef(sn1*rCellconv50b{dim1, sess, sc1}(trl, tpts), sn2*rCellconv50b{dim2, sess, sc2}(trl, tpts), 'rows', 'complete');
                            resCorr_wiTrialConv50b{dim1, dim2, sess, sc2}(trl) = R(1, 2);
                            
                            tpts2 = tpts(1 : 50 : end);
                            if length(tpts2) >= 4
                                [R, P] = corrcoef(sn1*rCellconv50b{dim1, sess, sc1}(trl, tpts2), sn2*rCellconv50b{dim2, sess, sc2}(trl, tpts2), 'rows', 'complete');
                                resCorr_wiTrialConv50bStep{dim1, dim2, sess, sc2}(trl) = R(1, 2);
                            end
                        end
                        
                    end
                    dataOut = r2z(resCorr_wiTrial{dim1, dim2, sess, sc2});
                    r = tanh(nanmean(dataOut));
                    meanCorrRorig{dim1, dim2, sess, sc2} = nanmean(resCorr_wiTrial{dim1, dim2, sess, sc2});
                    meanCorrR{dim1, dim2, sess, sc2} = r;
                    meanCorrZ{dim1, dim2, sess, sc2} = nanmean(dataOut);
                    
                    dataOut = r2z(resCorr_wiTrialConv{dim1, dim2, sess, sc2});
                    r = tanh(nanmean(dataOut));
                    meanCorrRorigC{dim1, dim2, sess, sc2} = nanmean(resCorr_wiTrialConv{dim1, dim2, sess, sc2});
                    meanCorrRconv{dim1, dim2, sess, sc2} = r;
                    meanCorrZconv{dim1, dim2, sess, sc2} = nanmean(dataOut);
                    
                    dataOut = r2z(resCorr_wiTrialConv50{dim1, dim2, sess, sc2});
                    r = tanh(nanmean(dataOut));
                    meanCorrRorigC50{dim1, dim2, sess, sc2} = nanmean(resCorr_wiTrialConv50{dim1, dim2, sess, sc2});
                    meanCorrRconv50{dim1, dim2, sess, sc2} = r;
                    meanCorrZconv50{dim1, dim2, sess, sc2} = nanmean(dataOut);
                    
                    dataOut = r2z(resCorr_wiTrialConv50b{dim1, dim2, sess, sc2});
                    r = tanh(nanmean(dataOut));
                    meanCorrRorigC50b{dim1, dim2, sess, sc2} = nanmean(resCorr_wiTrialConv50b{dim1, dim2, sess, sc2});
                    meanCorrRconv50b{dim1, dim2, sess, sc2} = r;
                    meanCorrZconv50b{dim1, dim2, sess, sc2} = nanmean(dataOut);
                    
                    dataOut = r2z(resCorr_wiTrialConv50bStep{dim1, dim2, sess, sc2});
                    r = tanh(nanmean(dataOut));
                    meanCorrRorigC50bStep{dim1, dim2, sess, sc2} = nanmean(resCorr_wiTrialConv50bStep{dim1, dim2, sess, sc2});
                    meanCorrRconv50bStep{dim1, dim2, sess, sc2} = r;
                    meanCorrZconv50bStep{dim1, dim2, sess, sc2} = nanmean(dataOut);
                end
            end
        end
    end
end


if doSave == 1
    save('E:\Matalb analyses\resCorr_wiTrial_230522', 'par_sess',...
        'resCorr_wiTrial', 'meanCorrRorig', 'meanCorrR', 'meanCorrZ', ...
        'filtW','resCorr_wiTrialConv', 'meanCorrRorigC', 'meanCorrRconv', 'meanCorrZconv',...
        'filtW50','resCorr_wiTrialConv50', 'meanCorrRorigC50', 'meanCorrRconv50', 'meanCorrZconv50',...
        'filtW50b','resCorr_wiTrialConv50b', 'meanCorrRorigC50b', 'meanCorrRconv50b', 'meanCorrZconv50b',...
        'resCorr_wiTrialConv50bStep', 'meanCorrRorigC50bStep', 'meanCorrRconv50bStep', 'meanCorrZconv50bStep');
%     save('E:\Matalb analyses\resCorr_wiTrial_230427', 'par_sess',...
%         'resCorr_wiTrial', 'meanCorrRorig', 'meanCorrR', 'meanCorrZ', ...
%         'filtW','resCorr_wiTrialConv', 'meanCorrRorigC', 'meanCorrRconv', 'meanCorrZconv',...
%         'filtW50','resCorr_wiTrialConv50', 'meanCorrRorigC50', 'meanCorrRconv50', 'meanCorrZconv50',...
%         'filtW50b','resCorr_wiTrialConv50b', 'meanCorrRorigC50b', 'meanCorrRconv50b', 'meanCorrZconv50b',...
%         'resCorr_wiTrialConv50bStep', 'meanCorrRorigC50bStep', 'meanCorrRconv50bStep', 'meanCorrZconv50bStep');
end






