% assignSignal
% Assign signal from sigSL and sigRL

clear sl rl ylims
switch dim
    case 'TinC'
        sl = sigSL.TinC*1000;
        rl = sigRL.TinC*1000;
        ylims = [0 30];
        if sub0 == 1
            ylims = [0 20];
        end
    case 'TinI'
        sl = sigSL.TinI*1000;
        rl = sigRL.TinI*1000;
        ylims = [0 20];
    case 'WhenD'
        sl = sigSL.whenC;
        rl = sigRL.whenC;
        ylims = [0 1.1];
    case 'PC1'
        sl = sigSL.PC1;
        rl = sigRL.PC1;
        ylims = [0 1.1];
        if sub0 == 1
            ylims = [0 0.7];
        end
    case 'MinC'
        sl = sigSL.MinC*1000;
        rl = sigRL.MinC*1000;
        ylims = [0 25];
    case 'MinI'
        sl = sigSL.MinI*1000;
        rl = sigRL.MinI*1000;
        ylims = [0 26];
    case 'WhatD'
        sl = sigSL.whatD;
        rl = sigRL.whatD;
        ylims = [-.1 1.1];
    case 'Ramp'
        sl = sigSL.ramp;
        rl = sigRL.ramp;
        ylims = [-.1 1.1];
        if sub0 == 1
            ylims = [0 0.6];
        end
    case 'WhenD_Tin'
        sl = sigSL.whenC_tin;
        rl = sigRL.whenC_tin;
        ylims = [0 1.1];
    case 'WhenD_noTin'
        sl = sigSL.whenC_noTin;
        rl = sigRL.whenC_noTin;
        ylims = [0 1.1];
    case 'ramp_noMin'
        sl = sigSL.ramp_noMin;
        rl = sigRL.ramp_noMin;
        ylims = [-0.2 1.1];
    case 'IMinReg'
        sl = sigSL.minReg *1000;
        rl = sigRL.minReg *1000;
        ylims = [-1 1];
        
end