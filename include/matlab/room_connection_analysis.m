%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paper: Structured Indoor Modeling
% published by: IEEE International Conference on Computer Vision (ICCV) 2015
% authors: Satoshi Ikehata, Hang Yan, Yasutaka Furukawa
% contact: sikehata@wustl.edu, yanhang@seas.wustl.edu,
% furukawa@wustl.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mergeList, doorList] = room_connection_analysis(WALLS, LABELS, PATH, Evidences, voxel_param, door_param, options)
% [mergeList, doorList] = ROOM_CONNECTION_ANALYSIS(WALLS, LABELS, PATH, Evidences, voxel_param, door_param, options)
% Analize the connection type of two wall nodes
% Merge:  then run room merging rule to merge two room nodes later
% Door: then run door addition rule to add doors later


display('[ROOM CONNECTION ANALYSIS]');
[WALLS2, LABELS2, PATH2, newClusterMap] = arrange_room_nodes(WALLS, LABELS, PATH, voxel_param);

propagated = geodesic_propagation(newClusterMap, double(Evidences.freespace_evidence>0));
[labels, connectivity] = find_connectivity(propagated);
if size(WALLS,1) ~= length(labels)
    mergeList = [];
    doorList = [];
    return;
end

doorList = check_wall_connection_type(connectivity, Evidences.freespace_evidence3d, Evidences.point_evidence3d, voxel_param, options,  door_param, WALLS2, LABELS2);

if isempty(cell2mat(doorList))
    doorList = {[]};
    mergeList = [];
    return;
end

% check door statistics

% door-> aisle
doorcount = 0;
doorsize = [];
indexdoor = [];
for i = 1 : length(doorList)
   if doorList{i}.aisleFlag == 0
       x0 = doorList{i}.x0;
       x1 = doorList{i}.x1;
       d = x1(1:2)-x0(1:2);
       d = d/norm(d);
       
       if d(1) == 1 | d(1) == -1                 
       doorsize(doorcount+1,1) = abs(x1(1) - x0(1));
       doorsize(doorcount+1,2) = x1(3) - x0(3);
       end
       
       if d(2) == 1 | d(2) == -1
       doorsize(doorcount+1,1) = abs(x1(2) - x0(2));
       doorsize(doorcount+1,2) = x1(3) - x0(3);
       end   
       indexdoor(doorcount+1,1) = i;
       doorcount = doorcount+1;

   end
end

xmax =min(doorsize(:,1)) + max(std(doorsize(:,1)) + 0.3*min(doorsize(:,1)));
for i = 1 : length(indexdoor)
   if doorsize(i,1) > xmax       
       doorList{indexdoor(i)}.aisleFlag = 1;
       display(sprintf('door %d is changed into door to aisle',i));
   end
end



% % aisle -> door
% 
% doorcount = 0;
% indexdoor = [];
% doorsize = [];
% 
% aislecount = 0;
% indexaisle = [];
% aislesize = [];
% 
% for i = 1 : length(doorList)
%    x0 = doorList{i}.x0;
%    x1 = doorList{i}.x1;  
%    d = x1(1:2)-x0(1:2);
%    d = d/norm(d);
%    if doorList{i}.aisleFlag == 0    
%        if d(1) == 1 | d(1) == -1                 
%        doorsize(doorcount+1,1) = abs(x1(1) - x0(1));
%        doorsize(doorcount+1,2) = x1(3) - x0(3);       
%        end       
%        if d(2) == 1 | d(2) == -1
%        doorsize(doorcount+1,1) = abs(x1(2) - x0(2));
%        doorsize(doorcount+1,2) = x1(3) - x0(3);
%        end      
%        indexdoor(doorcount+1,1) = i;
%        doorcount = doorcount + 1;
%    else % aisle
%        if d(1) == 1 | d(1) == -1                 
%        aislesize(aislecount+1,1) = abs(x1(1) - x0(1));
%        aislesize(aislecount+1,2) = x1(3) - x0(3);
%        end
%        
%        if d(2) == 1 | d(2) == -1
%        aislesize(aislecount+1,1) = abs(x1(2) - x0(2));
%        aislesize(aislecount+1,2) = x1(3) - x0(3);
%        end        
%        indexaisle(aislecount+1,1) = i;
%        aislecount = aislecount + 1;
%    end  
% end
% 
% for i = 1 : aislecount
%    for j = 1 : doorcount
%        if abs(aislesize(i,1)-doorsize(j,1))/options.leafsize(1) <= 2 && abs(aislesize(i,2)-doorsize(j,2))/options.leafsize(1) <= 2
% %            [abs(aislesize(i,1)-doorsize(j,1))/options.leafsize(1),  abs(aislesize(i,2)-doorsize(j,2))/options.leafsize(1)]
%            doorList{indexaisle(i)}.aisleFlag = 0;
%        end
%    end
% end

% figure, plot(doorsize(:,1), doorsize(:,2), '*'), axis equal
% axis([voxel_param.b0(1) voxel_param.b0(1)+voxel_param.bSize(1) voxel_param.b0(2) voxel_param.b0(2)+voxel_param.bSize(2)])


numRooms = length(WALLS2);
merge_matrix = zeros(numRooms, numRooms);
count_connectivity = zeros(numRooms, numRooms);
for i = 1 : size(doorList)
    if doorList{i}.aisleFlag == 1
        merge_matrix(min(doorList{i}.roomId1, doorList{i}.roomId2), max(doorList{i}.roomId1, doorList{i}.roomId2)) = 1;
    end    
end
[room1, room2] = ind2sub([numRooms, numRooms], find(merge_matrix==1));
mergeList = [room1 room2];
end



