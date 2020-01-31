function [POINT_room, CAM_room] = assign_room_point(POINT, CAM, CAMLIST, options, voxel_param, core_point_idx, mask)
% [POINT_room, CAM_room] = ASSIGN_ROOM_POINT(POINT, CAM, CAMLIST, options, voxel_param, core_point_idx, mask)
% Assign subset of the pointcloud to each room node based on the recovered room
% structure

indexList = setdiff(unique(mask),0);
%% assign points to core clusters
POINT_room = cell(1,1);
CAM_room = cell(1,1);

for i = 1 : length(indexList)
[P, C] = assign_point_room_labels(voxel_param, options, mask == indexList(i), POINT, CAM, core_point_idx);  
POINT_room{i,1} = P;
CAM_room{i,1} = C;
end
end