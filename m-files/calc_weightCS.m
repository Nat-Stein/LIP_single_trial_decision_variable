
% Not the latest version
load(fullfile(saveLoc, 'sig_allSessions'), 'sigSL')


ws{1} = sigSL.w.TinC;
ws{2} = sigSL.w.TinI;
ws{3} = sigSL.w.ramp;
ws{4} = sigSL.w.PC1;
ws{5} = sigSL.w.whatD;
ws{6} = sigSL.w.whenC;

allCs = nan(8,length(allW),length(allW));
for sess = 1 : 8
    for w1 = 1 : length(dNames)
        sc1 = 1; sn = 1; 
        for w2 = 1 : length(dNames)
            sc2 = 1; sn = 1; 
            x = []; y = [];
            if size(allW{w1}, 2) >= sess; x = allW{w1}{sc1, sess}; end
            if size(allW{w2}, 2) >= sess; y = allW{w2}{sc2, sess}; end
            if length(x) > 0 & length(y) > 0 & length(x) == length(y)
                xy   = dot(x,y);
                nx   = norm(x);
                ny   = norm(y);
                nxny = nx*ny;
                Cs   = xy/nxny;
                allCs(sess, w1, w2) = Cs;
            end
        end
    end
end


