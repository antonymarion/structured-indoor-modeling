function [WALLS2, LABELS2, PATH2] = delete_room_in_room(WALLS, LABELS, PATH, voxel_param)
% [WALLS2, LABELS2, PATH2] = DELETE_ROOM_IN_ROOM(WALLS, LABELS, PATH, voxel_param)
% Remove a room node if it is in the other room node coverage

indexList = [1:length(WALLS)]';
newClusterMap = path_to_mask(PATH, voxel_param);
indexValid = unique(newClusterMap);
WALLS2 = cell(1,1);
LABELS2 = cell(1,1);
PATH2 = cell(1,1);
count = 1;

for i = 1 : length(WALLS)
    if ~isempty(find(indexValid==i))
    WALLS2(count,1) = WALLS(i);
    LABELS2(count,1) = LABELS(i);
    PATH2(count,1) = PATH(i);
    count = count + 1;
    end
end


end

