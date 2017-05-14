clear;
%nx = 100;
length = 10;
step = 1;
window1 = 49;
window2 = 49+length;
while (window2 < 113)


%files{1} = './data/FIL_GoEFphi2_lfc.txt';
%files{2} = './data/FIL_GoEFqESV_lfc.txt';
%files{3} = './data/FIL_GoEFqI_lfc.txt';


%files{1} = './data/allPhi2_consistent_cleaned_LFC_avg_exp.txt'; %data for
%cell system paper
%files{2} = './data/allQESV_consistent_cleaned_LFC_avg_exp.txt';
%files{3} = './data/allQI_new_RAW3_adj_LFC_avg_consistent_NA_exp_PosONLY.txt';

files{1} = './data/IDMapping_consolidated_allPhi2_cleaned_lfc_avg.txt'; %semi-raw data
files{2} = './data/IDMapping_consolidated_allQESV_cleaned_LFC_avg.txt';
files{3} = './data/IDMapping_consolidated_allQI_new_RAW3_adj_LFC_avg.txt';

Par.numcluster = 10;   % number of clusters, needs to be specified.
Par.normalize  = 2;   % optional, 0 - no normalization; 1 - normalize to [0,1]; 2 - normalize to [0,0.5] and [0.5,1]. Default: 2
Par.anchor     = 100; %100; % optional, number of anchor points, default: 100
Par.maxit      = 500; % optional, maximum iterations of EM algorithm, default: 400
Par.Leps       = 0.1; %10;   % optional, termination criteria of EM algorithm, default: 1
Par.plot       = 0;   % optional, true - plot clustering result, default: true

Par.pca        = 0;   %draw the first two dimensions of pca or the original data, default: false
Par.start      = window1;   %49;%17; %1;   %start column position of each matrix
Par.end        = window2;  %112; %48; %16;  %end column position of each matrix
fiout=[int2str(window1) ,'.txt'];
Par.output     = ['./data/cluster',fiout]; % outout file name
disp(Par.output);
Idx = kde_em_clustering(files,Par);
window1 = window1 +step;
window2 = window2 +step;
end