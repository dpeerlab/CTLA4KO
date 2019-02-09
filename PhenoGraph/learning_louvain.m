a = cytof_data('~/Desktop/fcs_files/subsample_8k.fcs');
cyclinb1 = a.data(:, a.name_channel_map('cyclinb'));
ph3 = a.data(:, a.name_channel_map('ph3'));
prb = a.data(:, a.name_channel_map('prb'));
idu = a.data(:, a.name_channel_map('idu'));
data = [idu, cyclinb1, ph3, prb];
%%
[IDX,D] = knnsearch(data,data,'k',11,'distance','euclidean');
IDX(:, 1) = [];
D(:, 1) = [];
G = knn2jaccard(IDX);