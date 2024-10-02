function [R,ti] = bigSbigT(H,dt,tstep,twin,conditions,filter,min_num_bins)

% AZ
% function [R,ti] = bigSbigT(H,dt,tstep,twin,conditions,filter,min_msec)
% 
% Inputs
%     H: spike counts [trials x time steps]
%     dt: bin duration
%     tstep: time step for calculation of the rates
%     twin: time window for averaging
%     conditions: [trials x 1] or [trials x m]. Identifies the unique
%        conditions for averaging
%     filter: include/exclude trials
%     min_num_bins: minimum number of bins required to compute a rate;
%        otherwise, NaN

step = round(tstep / dt);
win = round(twin / dt);

[Y,~] = sum_steps_nan(H,step,win,2);
[N,ti] = sum_steps_nan(double(~isnan(H)),step,win,2);

if ~isvector(conditions)
    cond = nan(size(conditions,1),1);
    for i=1:size(conditions,2)
        cond(conditions(:,i)==1) = i;
    end
    conditions = cond;
end

uni_conditions = nanunique(conditions);

R = nan(length(uni_conditions),size(Y,2));
for i=1:length(uni_conditions)
    inds = conditions==uni_conditions(i) & filter==1;
    denom = nansum(N(inds,:));
    tind = denom >= min_num_bins;
    R(i,tind) = nansum(Y(inds,tind))./denom(tind);
end

% R = 1000/dt * R; %to firing rate

end