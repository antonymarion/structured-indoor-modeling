%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paper: Structured Indoor Modeling
% published by: IEEE International Conference on Computer Vision (ICCV) 2015
% authors: Satoshi Ikehata, Hang Yan, Yasutaka Furukawa
% contact: sikehata@wustl.edu, yanhang@seas.wustl.edu,
% furukawa@wustl.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [WALLS2, LABELS2, PATH2, doorList2, mergeList2] = delete_empty_rooms(WALLS, LABELS, PATH, voxel_param, doorList, unexplained_indeces)
% [WALLS2, LABELS2, PATH2, doorList2, mergeList2] = DELETE_EMPTY_ROOMS(WALLS, LABELS, PATH, voxel_param, doorList, unexplained_indeces)
% Remove room nodes if they have no connection to other room nodes

if isempty(unexplained_indeces)
    
WALLS2 = cell(1,1);
LABELS2 = cell(1,1);
PATH2 = cell(1,1);
doorList2 = doorList;
if isempty(doorList)
    WALLS2 = WALLS;
    LABELS2 = LABELS;
    PATH2 = PATH;
    return;
end
numRooms = length(WALLS);
connectivity_index = zeros(numRooms, numRooms);
for i = 1 : length(doorList)
    if doorList{i}.aisleFlag == 0
    room1 = doorList{i}.roomId1;
    room2 = doorList{i}.roomId2;
    connectivity_index(room1,room2) = 1;
    connectivity_index(room2,room1) = 1;    
    end
end
 
valid_indices = find(sum(connectivity_index)>0);
invalid_indices = find(sum(connectivity_index)==0);
% update roomid
before = [1:length(WALLS)]';
after = zeros(length(WALLS),1);
count = 1;
for i = 1 : length(valid_indices)
    after(valid_indices(i)) = count;
    WALLS2{count,1} = WALLS{valid_indices(i)};
    LABELS2{count,1} = LABELS{valid_indices(i)};
    PATH2{count,1} = PATH{valid_indices(i)};
    count = count + 1;
end

for i = 1 : length(doorList)
    doorList2{i}.roomId1 = after(doorList{i}.roomId1);
    doorList2{i}.roomId2 = after(doorList{i}.roomId2);
end

else
WALLS2 = cell(1,1);
LABELS2 = cell(1,1);
PATH2 = cell(1,1);
doorList2 = cell(1,1);
if isempty(doorList)
    WALLS2 = WALLS;
    LABELS2 = LABELS;
    PATH2 = PATH;
    return;
end
numRooms = length(WALLS);
connectivity_index = zeros(numRooms, numRooms);
for i = 1 : length(doorList)
    room1 = doorList{i}.roomId1;
    room2 = doorList{i}.roomId2;
    connectivity_index(room1,room2) = 1;
    connectivity_index(room2,room1) = 1;    
end
 
invalid_indices = find(sum(connectivity_index)<=1);
invalid_indices = intersect(unexplained_indeces, invalid_indices);
valid_indices = setdiff([1:length(WALLS)]', invalid_indices);

count = 0;
for i = 1: length(doorList)
    if isempty(find(invalid_indices==doorList{i}.roomId1)) && isempty(find(invalid_indices==doorList{i}.roomId2))
        doorList2{count+1,1} = doorList{i};
        count = count + 1;
    end  
end

% update roomid
before = [1:length(WALLS)]';
after = zeros(length(WALLS),1);
count = 1;
for i = 1 : length(valid_indices)
    after(valid_indices(i)) = count;
    WALLS2{count,1} = WALLS{valid_indices(i)};
    LABELS2{count,1} = LABELS{valid_indices(i)};
    PATH2{count,1} = PATH{valid_indices(i)};
    count = count + 1;
end

for i = 1 : length(doorList2)
    doorList2{i,1}.roomId1 = after(doorList2{i,1}.roomId1);
    doorList2{i,1}.roomId2 = after(doorList2{i,1}.roomId2);
end

numRooms = length(WALLS2);
merge_matrix = zeros(numRooms, numRooms);

for i = 1 : size(doorList2)
    if doorList2{i}.aisleFlag == 1
        merge_matrix(min(doorList2{i}.roomId1, doorList2{i}.roomId2), max(doorList2{i}.roomId1, doorList2{i}.roomId2)) = 1;
    end    
end
[room1, room2] = ind2sub([numRooms, numRooms], find(merge_matrix==1));
mergeList2 = [room1 room2];
end

end

