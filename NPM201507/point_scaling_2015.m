function D_scale = point_scaling_2015(D,type)

num_plant = length(D);
points = [];
min_pos = 0;%0.25;
max_neg = 0; %-0.25;

for i = 1:num_plant
    points = [points, D{i}];
end

if type == 1
    [~,PS1] = mapminmax(points,0,1);
    for i = 1:num_plant
        D_scale{i} = mapminmax('apply',D{i},PS1);
    end
elseif type == 2
    points_pos = points;
    points_neg = points;
    points_pos(points_pos<0) = 0;
    points_neg(points_neg>0) = 0;
    [~,PS1] = mapminmax(points_pos,min_pos,0.5);
    [~,PS2] = mapminmax(points_neg,-0.5,max_neg);
    
    for i = 1:num_plant
        D_pos = D{i};
        D_neg = D{i};
        D_pos(D_pos<0) = 0;
        D_neg(D_neg>0) = 0;
        D_scale{i} = mapminmax('apply',D_pos,PS1) + mapminmax('apply',D_neg,PS2); %D_pos,PS1);
    end
    
end