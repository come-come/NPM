function [Idx] = kde_em_clustering(files, Par)

%set default values for the parameters that are not specified by user
if ~isfield(Par,'normalize')
    Par.normalize  = 2;
end
if ~isfield(Par,'plot')
    Par.plot  = 1;
end

%load data from files %assume the row order are the same for all input
%files
for m = 1:length(files)
    FF = importdata(files{m});
    FF.data = FF.data(:, Par.start:Par.end);
    [num_plant,num_timepoint] = size(FF.data);
    num_plant = num_plant - 1;
    name_plant = FF.rowheaders(2:num_plant+1);
    data_meas{m} = FF.data(2:num_plant+1,:);
end

names = FF.textdata;

%delete elements with NAN
for i = 1:num_plant
    datalines = [];
    for m = 1:length(files)
        datalines = [datalines;data_meas{m}(i,:)];
    end
    D{i} = datalines;
    [row, col] = find(isnan(D{i}));  % delete columns with NAN
    D{i}(:,col) = [];
%    [row, col] = find(D{i}==0); % delete columns with 0
%    D{i}(:,col) = [];
end

%data normalization/scaling
D_origin = D;
if Par.normalize
    D = point_scaling_2015(D,Par.normalize);
end

%clustering on the normalized data
kkk = Par.numcluster;
for i = 1:10
    Idx = kde_tm_EM_2015(D,kkk,Par);
    if length(unique(Idx))==kkk
        break;
    end
    fprintf('Cannot find enough clusters. Trying again.');
end

if length(unique(Idx)) < kkk
    fprintf('Cannot find enough clusters. Please reduce number of clusters.');
    return
end

%dlmwrite(Par.output,Idx); %only output cluster IDs

fileID = fopen(Par.output,'w'); %output both row ID and the corresponding cluster ID
formatSpec = '%s	%d	%d\n';
NewIdx = cat(1,1,Idx);
[nrows,~] = size(NewIdx);
for row = 1:nrows
    fprintf(fileID,formatSpec,names{row}, NewIdx(row), 2);
end
fclose(fileID);

if Par.plot
    if Par.pca 
        if size(D_origin{1},1) > 2
            D_pca = point_pca_2015(D_origin,2);
        else
            D_pca = D_origin;
        end
    else
        D_pca = D_origin;
    end
    
    plot_clusters_2015(D_pca,Idx,kkk);
    
end
