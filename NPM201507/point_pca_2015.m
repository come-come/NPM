function D_pca = point_pca_2015(D,dim)

num_plant = length(D);
points = [];
for i = 1:num_plant
    points = [points, D{i}];
end

[U,S,V] = svd(cov(points'));
points_pca = (points'*U(:,1:dim))';

cum = 1;
for i = 1:num_plant
    D_pca{i} = points_pca(:,cum:cum+size(D{i},2)-1);
    cum = cum+size(D{i},2);
end