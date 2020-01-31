function [Pout, Cout] = assign_point_room_labels(voxel_param, options, mask, POINT, CAM, core_point_idx)
% [Pout, Cout] = ASSIGN_POINT_ROOM_LABELS(voxel_param, options, mask, POINT, CAM, core_point_idx);
numPoints = size(POINT,1);
P = cell2mat(POINT);
C = cell2mat(CAM);

if ~isempty(core_point_idx)
P = P(core_point_idx,:);
C = C(core_point_idx);
end

X = floor((P(:,3)-voxel_param.b0(1))/options.leafsize(1))+1;
Y = floor((P(:,4)-voxel_param.b0(2))/options.leafsize(2))+1;
validindex = find(X>=1 & X<=voxel_param.gNum(1) & Y>=1 & Y<= voxel_param.gNum(2));

ID = (X-1)*voxel_param.gNum(2) + Y;
validID = ID(validindex);
indexRoom = find(mask(validID)>0);
Pout = P(validindex(indexRoom),:);
Cout = C(validindex(indexRoom),:);

end

