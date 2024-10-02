% general_settings for LIP single trial paper



xLimsSL = [0 0.6];
xLimsRL = [-0.3 0.1];
xLimsView = [0 0.4];

fontSize = 6;
sFont=6;


colorsYeBl = {{[102 51 0]/255, [204 102 0]/255, [255 128 0]/255,...
    [255 178 102]/255, [255 229 204]/255, 0.6118*[1 1 1]},...
    {[0 51 102]/255, [0 102 204]/255, [0 128 255]/255,...
    [102 178 255]/255, [153 204 255]/255, 0.6118*[1 1 1]}};
clear colorsYeBl_ord
for i = 1 : 6
    colorsYeBl_ord{i} = colorsYeBl{2}{i};
    colorsYeBl_ord{12-i} = colorsYeBl{1}{i};
end



cutoff = 0.7; ls = '-'; lw = 1;
smooth_win_t = 0.025;

fw = 50; filtW50b = ones(1, fw)/fw; % used for within trial correlations

