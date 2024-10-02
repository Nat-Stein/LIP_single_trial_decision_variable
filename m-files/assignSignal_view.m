% assignSignal
% Assign signal from sigSL and sigRL


clear sl rl ylims
switch dim
    case 'TinC'
        sl = sigSLv.TinC*1000;
        ylims = [0 30];
        if sub0 == 1
            ylims = [0 20];
        end
    case 'TinI'
        sl = sigSLv.TinI*1000;
        ylims = [0 20];
    case 'WhenD'
        sl = sigSLv.whenC;
        ylims = [0 1.1];
    case 'PC1'
        sl = sigSLv.PC1;
        ylims = [0 1.1];
        if sub0 == 1
            ylims = [0 0.6];
        end
    case 'MinC'
        sl = sigSLv.MinC*1000;
        ylims = [0 25];
    case 'MinI'
        sl = sigSLv.MinI*1000;
        ylims = [0 30];
    case 'WhatD'
        sl = sigSLv.whatD;
        ylims = [-.1 1.1];
    case 'Ramp'
        sl = sigSLv.ramp;
        ylims = [-.1 1.1];
        if sub0 == 1
            ylims = [0 0.6];
        end
    case 'WhenD_Tin'
        sl = sigSLv.whenC_tin;
        ylims = [0 1.1];
    case 'WhenD_noTin'
        sl = sigSLv.whenC_noTin;
        ylims = [0 1.1];
end