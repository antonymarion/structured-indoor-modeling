%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paper: Structured Indoor Modeling
% published by: IEEE International Conference on Computer Vision (ICCV) 2015
% authors: Satoshi Ikehata, Hang Yan, Yasutaka Furukawa
% contact: sikehata@wustl.edu, yanhang@seas.wustl.edu,
% furukawa@wustl.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [walls, label, sPath] = room_wall_path_extraction_normal(xStart, xGoal, point_evidence, freespace_evidence, wall_param, voxel_param, options, coremask, nx, ny)
% [walls, label, sPath] = ROOM_WALL_PATH_EXTRACTION_NORMAL(xStart, xGoal, point_evidence, freespace_evidence, wall_param, voxel_param, options, coremask, nx, ny)
% Solve modified piecewise planar reconstruction using dynamic programming
% Piecewise Planar and Compact Floorplan Reconstruction from Images Ricardo
% Cabral and Yasutaka Furukawa CVPR 2014.

% perform filtering 
kernel = fspecial('gaussian', 2*ceil(wall_param.margin_size/options.leafsize(1))+1, wall_param.margin_sigma);

point_evidence = imfilter(point_evidence, kernel);
point_evidence = normalize_evidence(point_evidence, -1);
freespace_evidence = normalize_evidence(freespace_evidence, -1);

% run shortest path reconstruction
[sPath, min_cost, path_length, count] = mex_shortest_path_reconstruction_normal(point_evidence, coremask, freespace_evidence, xStart, xGoal, wall_param.thresh, wall_param.penalty, 0, wall_param.margin_sigma, wall_param.min_path_length/options.leafsize(1), wall_param.freespace_weight, nx(:), ny(:));

% upgrage small segment
if length(sPath) < 4
    walls = [];
    label = [];
    sPath = []; 
    warning('Can not find the path')
    return;
end


[y x] = ind2sub([voxel_param.gNum(2), voxel_param.gNum(1)],sPath+1);

% remove start end 
valid_ids = ones(length(sPath), 1);
if (x(end) == x(1) & x(1) == x(2)) |  (y(end) == y(1) & y(1) == y(2))
    valid_ids(1) = 0;
end

if (x(end-1) == x(end) & x(end) == x(1)) |  (y(end-1) == y(end) & y(end) == y(1))
    valid_ids(end) = 0;
end

for i = 1 : length(sPath)-2
    if (x(i) == x(i+1) & x(i+1) == x(i+2)) |  (y(i) == y(i+1) & y(i+1) == y(i+2))
        valid_ids(i+1) = 0;
    end        
end

sPath = sPath(valid_ids==1);
[y x] = ind2sub([voxel_param.gNum(2), voxel_param.gNum(1)],sPath+1);
xx = voxel_param.bSize(1)*(x-1)/voxel_param.gNum(1)+voxel_param.b0(1);
yy = voxel_param.bSize(2)*(y-1)/voxel_param.gNum(2)+voxel_param.b0(2);
walls = [xx yy];
label = zeros(size(walls,1),1);
xx = x;
yy = y;
xx(end+1) = xx(1);
yy(end+1) = yy(1);
for i = 1 : size(walls,1)
    if xx(i+1)-xx(i) > 0 
        label(i) = 1;
    end
    
    if xx(i+1)-xx(i) < 0 
        label(i) = -1;
    end
    
    if yy(i+1)-yy(i) > 0 
        label(i) = 2;
    end
    
    if yy(i+1)-yy(i) < 0 
        label(i) = -2;
    end
end



end
