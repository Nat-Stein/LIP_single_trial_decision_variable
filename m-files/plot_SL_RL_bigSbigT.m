function [p,lineas,out] = plot_SL_RL_bigSbigT(s,smooth_win_t,conditions,filter,doPlot,...
    legend_str, varargin)
% function [p,Ton,Son,Toff,Soff,filter] = plot_lock_on_off(s,sm,conditions,filter,doPlot)

% INPUTS
% s is a struct: 
%         s.g_on: spike counts, aligned on dot onset [trials x time steps]
%         s.g_off: spike counts, aligned on RT [trials x time steps]
%         s.t_on: time vector for s.g_on
%         s.t_off: time vector for s.g_off
% smooth_win_t: time window for computing the rates
% conditions: either a vector identifying the unique conditions, or a
% matrix of 1s and 0s in which each column is a different condition [trials
% x 1, or trials x conditions]
% filter: vector of 1s and 0s indiciating which trials to include in the
% analysis [trials x 1]
% doPlot: boolean, to plot or not 
% legend_str: a cell with the name of each unique condition, to add a
% legend to the plot (optional; can be an empty cell) 


cutoff = 0.7;
independent_cutoff = false;
p_handle = [];
colortype = 0;
for i=1:2:length(varargin)
    if isequal(varargin{i},'cutoff')
        cutoff = varargin{i+1};
    elseif isequal(varargin{i},'independent_cutoff')
        independent_cutoff = varargin{i+1};
    elseif isequal(varargin{i},'p_handle')
        p_handle = varargin{i+1};
    elseif isequal(varargin{i},'color')
        colortype = varargin{i+1};
    elseif isequal(varargin{i},'colors')
        colores = varargin{i+1};
    elseif isequal(varargin{i},'min_num_bins')
        min_num_bins = varargin{i+1};
    end
end

if nargin<3 || isempty(filter)
    filter = ones(size(conditions,1),1)==1;
end

if nargin<5 || isempty(doPlot)
    doPlot = 1;
end
if nargin<6
    legend_str = [];
end


lineas = [];
out = [];

if isvector(conditions)
    conditions = adummyvar(conditions);
end

ncond = size(conditions,2);
if ~exist('colores','var')
    if ncond>1
        colores = rainbow_colors(ncond,'colortype',colortype);
        % colores = movshon_colors(ncond);
    else
        colores = [0 0 1];
    end
end


%ignore from cutoff rows that are all NaN
mn = mean(all(isnan(s.g_on')),2);
cutoff = mn + (1-mn)*cutoff;

if (doPlot>0)
    if ~isempty(p_handle)
        p = p_handle;
    else
        p = publish_plot(1,2);
    end
else
    p = [];
end

% bigSbigT
dt = s.t_on(2)-s.t_on(1);
% min_num_bins = 50; % read from elsewhere...
[Son,ti_on] = bigSbigT(s.g_on,dt,dt,smooth_win_t,conditions,filter,min_num_bins);
ton = s.t_on(ti_on);


dt = s.t_off(2)-s.t_off(1);
[Soff,ti_off] = bigSbigT(s.g_off,dt,dt,smooth_win_t,conditions,filter,min_num_bins);
toff = s.t_off(ti_off);


for i = 1:ncond

    son = Son(i,:);
    soff = Soff(i,:);

    ind = conditions(:,i)==1 & filter==1;
    
%     [~,son,erron ] = curva_media(Son, ones(size(filter)) ,ind,0);
%     [~,soff,erroff] = curva_media(Soff, ones(size(filter)),ind,0);


    if (independent_cutoff)
        nans_on  = nanmean(isnan(s.g_on(ind,ti_on)));
        nans_off = nanmean(isnan(s.g_off(ind,ti_off)));
    else
        nans_on  = nanmean(isnan(s.g_on(:,ti_on)));
        nans_off = nanmean(isnan(s.g_off(:,ti_off)));
    end

    tind_on  = nans_on<cutoff;
    tind_off = nans_off<cutoff;


    I1 = tind_on & ~isnan(son);
    out(i).on.t = ton(I1);
    out(i).on.y = son(I1);

    I2 = tind_off & ~isnan(soff);
    out(i).off.t = toff(I2);
    out(i).off.y = soff(I2);


    if (doPlot>0)
        p.next();

        plot(out(i).on.t, out(i).on.y,'color',colores(i,:));
        hold all

        p.next();

        lineas(i) = plot(out(i).off.t, out(i).off.y,'color',colores(i,:));
        hold all

    end
end

if (doPlot>0)
    same_ylim(p.h_ax);

    if isempty(legend_str)
        % hl = legend_n(1:ncond,'hline',lineas);
        hl = [];
    else
        hl = legend(lineas,legend_str);
    end
    if ~isempty(hl)
        set(hl,'location','best','interpreter','none')
    end

    p.format('FontSize',15)

    %     set(hl,'FontSize',8)

    set(p.h_ax(2),'YAxisLocation','right')

    %ha = brokenAxis([0.5 0.11 0.012 0.06],'vertical');
    %set(gcf,'CurrentAxes',p.h_ax(2));

    set(gcf,'Position',[457  399  872  333])
end

end




function [R,ti] = bigSbigT(H,dt,tstep,twin,conditions,filter,min_num_bins)

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