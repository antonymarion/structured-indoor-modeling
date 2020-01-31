%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paper: Structured Indoor Modeling
% published by: IEEE International Conference on Computer Vision (ICCV) 2015
% authors: Satoshi Ikehata, Hang Yan, Yasutaka Furukawa
% contact: sikehata@wustl.edu, yanhang@seas.wustl.edu,
% furukawa@wustl.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [WALLS2, LABELS2, PATH2] = merge_rooms(mergeList, doorList, WALLS, LABELS, PATH, Evidences, wall_param, voxel_param, options)
% [WALLS2, LABELS2, PATH2] = MERGE_ROOMS(mergeList, doorList, WALLS, LABELS, PATH, Evidences, wall_param, voxel_param, options)
% Apply room merging rule on two room nodes that are classified as MERGE connection at the previous stage 
% After the room merging, we also apply room reconstruction rule on the merged room node 


cmask = Evidences.freespace_evidence > 0;
cmask = imopen(cmask, strel('disk',1)); % for robustness (maybe cannot work for sparse point evidence...)


display('[ROOM MERGING]')
display(sprintf('We found %d merge connections...', size(mergeList,1)));
newClusterMap = path_to_mask(PATH, voxel_param);
% propagated = geodesic_propagation(newClusterMap, Evidences.binarymask);
propagated = geodesic_propagation(newClusterMap, cmask);
indexList = [1:length(unique(newClusterMap))-1];

isValid = zeros(length(indexList), 1);
for i = 1 : length(doorList)
    if doorList{i}.aisleFlag == 1
    isValid(doorList{i}.roomId1) = 1;
    isValid(doorList{i}.roomId2) = 1;
    end
end

if isempty(mergeList)
   newClusterMap2 =  newClusterMap;
   propagated2 = propagated;
   return;
end
newClusterMap2 = newClusterMap;
propagated2 = propagated;

core_indeces = find(isValid==1);

%% computing the coremask

mergemask = zeros(voxel_param.gNum(2),voxel_param.gNum(1));
coremask = zeros(voxel_param.gNum(2),voxel_param.gNum(1));
display('Computing the core free-space evidence...')
for i = 1 : length(core_indeces)
    index_cluster = core_indeces(i);
    supportmask = propagated == index_cluster;
    supportmask(find(cmask==0)) = 0;
    inv_mask = 1-supportmask;    
    [yy xx] = ind2sub([voxel_param.gNum(2), voxel_param.gNum(1)], find(supportmask>0));
    integral_inv_mask = integralImage(inv_mask);
    integral_mask_pt = integralImage(Evidences.point_evidence_wall);
    X = mex_compute_core_freespace(double(supportmask), double(integral_inv_mask), double(integral_mask_pt), min(xx)-1, max(xx)-1, min(yy)-1, max(yy)-1); % matlab coordinate
    [XX YY] = meshgrid([X(1):X(3)],  [X(2):X(4)]);    
    index = (XX-1)*voxel_param.gNum(2) + YY;
    mergemask(index) = index_cluster;      
    
    [XX2 YY2] = meshgrid([X(1)+wall_param.core_freespace_margin:X(3)- wall_param.core_freespace_margin],  [X(2)+wall_param.core_freespace_margin:X(4)- wall_param.core_freespace_margin]);    
    index2 = (XX2-1)*voxel_param.gNum(2) + YY2;
    coremask(index2) = index_cluster;          
    display(sprintf('...Size of CFE (room %d) = %d', i, length(index)));
end



%% delete empty core cluster
emptyClusterList = setdiff(core_indeces,unique(coremask));
indexList(emptyClusterList) = -1;  

badIndex = [];
while ~isempty(mergeList)    
    area1 = zeros(length(mergeList),1);
    area2 = zeros(length(mergeList),1);
        
    for i = 1 : size(mergeList,1)
        area1(i) = length(find(newClusterMap2 == indexList(mergeList(i, 1))));
        area2(i) = length(find(newClusterMap2 == indexList(mergeList(i, 2))));       
    end      
    
    [~, idx] = max(area1+area2);
    index1 = mergeList(idx,1);
    index2 = mergeList(idx,2);
    mergeList = mergeList(setdiff([1:size(mergeList,1)], idx),:);
    
    if indexList(index1) == indexList(index2)  | indexList(index1) == -1 | indexList(index2) == -1 % continue if the core mask is empty
        continue;
    end
    display(sprintf('Merge core regions (%d, %d)',indexList(index1),indexList(index2)));

    for i = 1 : length(doorList)
       if doorList{i}.roomId1 == index1 & doorList{i}.roomId2 == index2
           doorIdx = i;
       end
    end
      
    newClusterMap2(newClusterMap2==indexList(doorList{doorIdx}.roomId1)) = indexList(index1);
    newClusterMap2(newClusterMap2==indexList(doorList{doorIdx}.roomId2)) = indexList(index1);  
    propagated2(propagated2==indexList(doorList{doorIdx}.roomId1)) = indexList(index1);
    propagated2(propagated2==indexList(doorList{doorIdx}.roomId2)) = indexList(index1);   
    mask1 = coremask == indexList(doorList{doorIdx}.roomId1);
    mask2 = coremask == indexList(doorList{doorIdx}.roomId2);
    fs = propagated2 == indexList(index1);    
    
    walls1 = WALLS{doorList{doorIdx}.roomId1};
    walls1(end+1,:) = walls1(1,:);
    x1s = walls1(doorList{doorIdx}.wallId1,:);
    x1e = walls1(doorList{doorIdx}.wallId1+1,:);
    
    walls2 = WALLS{doorList{doorIdx}.roomId2};
    walls2(end+1,:) = walls2(1,:);
    x2s = walls2(doorList{doorIdx}.wallId2,:);
    x2e = walls2(doorList{doorIdx}.wallId2+1,:);
    x = floor(voxel_param.gNum(1)*([x1s(1);x1e(1);x2s(1);x2e(1)]-voxel_param.b0(1))/voxel_param.bSize(1))+1;
    y = floor(voxel_param.gNum(2)*([x1s(2);x1e(2);x2s(2);x2e(2)]-voxel_param.b0(2))/voxel_param.bSize(2))+1;
    x(end+1) = x(1);
    y(end+1) = y(1);    
% 
    mask = mex_fill_inner_rectangle(x, y, voxel_param.gNum(1), voxel_param.gNum(2)); % input is matlab coordinate
    mask = reshape(mask, voxel_param.gNum(2), voxel_param.gNum(1));
    
    mask(newClusterMap==doorList{doorIdx}.roomId1) = 1;
    mask(newClusterMap==doorList{doorIdx}.roomId2) = 1;   
    
    boundary_r1 = mex_findBoundary(voxel_param.gNum(1), voxel_param.gNum(2), double(mask1(:)));
    boundary_indeces_r1 = find(boundary_r1 == 1);
    [yBoundary_r1 xBoundary_r1] = ind2sub([voxel_param.gNum(2), voxel_param.gNum(1)],boundary_indeces_r1);

    boundary_r2 = mex_findBoundary(voxel_param.gNum(1), voxel_param.gNum(2), double(mask2(:)));
    boundary_indeces_r2 = find(boundary_r2 == 1);
    [yBoundary_r2 xBoundary_r2] = ind2sub([voxel_param.gNum(2), voxel_param.gNum(1)],boundary_indeces_r2);

    sList = [xBoundary_r1, yBoundary_r1]-1;
    
    invalidmask = (cmask == 0);
    invalidmask = imdilate(invalidmask,strel('disk',10));
    

        mask(invalidmask == 1) = 0;        
        
        flag = -1;
        CC = bwconncomp(mask, 8);
        for i = 1 : CC.NumObjects
            test = zeros(size(mask));
            test(CC.PixelIdxList{i}) = 1;
            if ~isempty(intersect((sList(1:5:end,1)-1)*voxel_param.gNum(2)+sList(1:5:end,2), find(test==1))) & ~isempty(intersect(find(mask2==1), find(test==1)))
                flag = 1; 
                break;
            end
        end

%         figure, imagesc(mergemask), axis equal
        if flag == 1
        sPath = mex_shortest_path(sList(1:5:end,:)', double(mask2), double(mask));
        mergemask(mergemask==indexList(doorList{doorIdx}.roomId1)) = indexList(index1);
        mergemask(mergemask==indexList(doorList{doorIdx}.roomId2)) = indexList(index1);
        coremask(coremask==indexList(doorList{doorIdx}.roomId1)) = indexList(index1);
        coremask(coremask==indexList(doorList{doorIdx}.roomId2)) = indexList(index1);
        coremask(sPath==1) = indexList(index1);    
        indexList(find(indexList == indexList(index2))) = indexList(index1);
        else
        warning('cannot find the good connection')
        end  
end

% figure('WindowStyle', 'docked', 'Name', 'mergemask'), imagesc(mergemask), axis equal;
% figure('WindowStyle', 'docked', 'Name', 'coremask'), imagesc(coremask), axis equal;
% figure('WindowStyle', 'docked', 'Name', 'propagated'), imagesc(propagated2), axis equal

%% update path

indexall = setdiff(unique(mergemask),-1);
for i = 1 : length(indexall)
   if length(find(indexall(i)==indexList))==1
      emptyClusterList(end+1) = indexall(i);
   end
end

core_indeces = [];
unique_indeces = setdiff(unique(indexList),-1);
for i = 1 : length(emptyClusterList) 
    mergemask(mergemask==emptyClusterList(i)) = 0; 
    propagated2(propagated2==emptyClusterList(i)) = 0;
    coremask(coremask==emptyClusterList(i)) = 0;
end


count = 1;
for i = 1 : length(unique_indeces)
   if length(find(indexList == unique_indeces(i))) >= 2
       core_indeces(count) = unique_indeces(i);
       count = count + 1;
   end
end

index = setdiff(unique(coremask),0);  

numcluster = length(core_indeces);
num_connectivity = zeros(numcluster, 1);
[labels connectivity]  = find_connectivity(propagated2);
for i = 1 : length(core_indeces)
num_connectivity(i) = length( connectivity{find(labels == core_indeces(i))} );
end
order_clusters = sortrows([[1:numcluster]',num_connectivity], 2);

constraintmask = zeros(voxel_param.gNum(2), voxel_param.gNum(1));
count = 1;

for i = 1 : length(unique_indeces)
    if isempty(intersect(unique_indeces(i), core_indeces)) & isempty(intersect(unique_indeces(i), find(isValid==1)))
    constraintmask(find(newClusterMap2==unique_indeces(i))) = 1;  
    WALLS2{count,1} = WALLS{unique_indeces(i)};
    LABELS2{count,1} =  LABELS{unique_indeces(i)};
    PATH2{count,1} =  PATH{unique_indeces(i)};
    count = count+1;
    end
end

propagated2 = zeros(size(propagated));
for i = 1 : length(indexList)
    propagated2(find(propagated==i)) = indexList(i);
end

for i = 1 : length(core_indeces)
display(sprintf('Room Cluster %d: Shortest path reconstruction (# room is %d)', i, length(core_indeces)));   
index_cluster = core_indeces(order_clusters(i,1));
supportmask = propagated2 == index_cluster;
invalidmask = double(mergemask ~= index_cluster & mergemask> 0);    

constraintmask(cmask==0) = 1;


point_evidence_cluster = zeros(size(Evidences.point_evidence_wall));
point_evidence_cluster(propagated2==index_cluster) = Evidences.point_evidence_wall(propagated2 == index_cluster);

point_mask = mergemask == index_cluster;
boundary_point = mex_findBoundary(voxel_param.gNum(1), voxel_param.gNum(2), double(point_mask(:)));
boundary_indeces_point = find(boundary_point == 1);
[yBoundary_point xBoundary_point] = ind2sub([voxel_param.gNum(2) voxel_param.gNum(1)],boundary_indeces_point);

core_room = double(coremask==index_cluster);    
points = [xBoundary_point(:), yBoundary_point(:)];
path_evidence = zeros(size(points,1),1);

keepoutmask = zeros(size(cmask));
keepoutmask(constraintmask==1) = 1;
keepoutmask(invalidmask==1) = 1;
keepoutmask = imdilate(keepoutmask, strel('disk', 1));
[xStart, xGoal, core_room] = mex_find_start_end(core_room, points'-1,  double(path_evidence), double(point_evidence_cluster), double(keepoutmask(:))); 
core_room(constraintmask==1) = 1;
core_room(invalidmask==1) = 1;   


display('... Running 1st Dijkstra algorithm')
% tic;
[walls, label, sPath] = room_wall_path_extraction_normal(xStart, xGoal, point_evidence_cluster, Evidences.freespace_evidence_sparse, wall_param, voxel_param, options, core_room, Evidences.nx, Evidences.ny);   
% time = toc;
% display(sprintf('Elapsed time is %.02f sec', toc))


%% path evaluation and find new start-end point
[y x] = ind2sub([voxel_param.gNum(2), voxel_param.gNum(1)], [sPath(:)+1;sPath(1)+1]);    
points = cell(length(sPath), 1);
path_evidence = cell(length(sPath), 1);
for i = 1 : length(sPath)
    d = [x(i+1)-x(i),y(i+1)-y(i)];
    numpix = norm(d)+1;
    d = d/(numpix-1);
    X = [x(i)+([1:numpix]-1)*d(1);y(i)+([1:numpix]-1)*d(2)]';
    index_path = (X(:,1)-1)*voxel_param.gNum(2) + X(:,2);
    [yy xx] = ind2sub([voxel_param.gNum(2), voxel_param.gNum(1)], index_path);
    points{i,1} = [xx,yy];
    path_evidence{i,1} = sum(Evidences.point_evidence_wall(index_path))*ones(length(xx),1);
end

points = cell2mat(points);
path_evidence = cell2mat(path_evidence);
core_room = double(coremask==index_cluster);    

keepoutmask = zeros(size(core_room));
keepoutmask(constraintmask==1) = 1;
keepoutmask(invalidmask==1) = 1;
keepoutmask = imdilate(keepoutmask, strel('disk', 1));
[xStart, xGoal, core_room] = mex_find_start_end(core_room, points'-1,  double(path_evidence), double(point_evidence_cluster), double(keepoutmask)); 
core_room(constraintmask==1) = 1;
core_room(invalidmask==1) = 1;   


display('... Running 2nd Dijkstra algorithm')
% tic;
[walls, label, sPath] = room_wall_path_extraction_normal(xStart, xGoal, point_evidence_cluster, Evidences.freespace_evidence_sparse, wall_param, voxel_param, options, core_room, Evidences.nx, Evidences.ny);   
% time =  toc;
% display(sprintf('Elapsed time is %.02f sec', toc))

kernel = fspecial('gaussian', 2*ceil(wall_param.margin_size/options.leafsize(1))+1, wall_param.margin_sigma);
test = imfilter(point_evidence_cluster, kernel);
test = normalize_evidence(test, -1);
figure('WindowStyle', 'docked', 'Name', sprintf('Room %d',i)), imagesc(reshape(test, voxel_param.gNum(2), voxel_param.gNum(1))), axis equal
hold on
plot(xStart(1)+1,xStart(2)+1,'or','MarkerEdgeColor','g','MarkerSize',6)
plot(xGoal(1)+1,xGoal(2)+1,'*r','MarkerEdgeColor','b','MarkerSize',6)
[y x] = ind2sub([voxel_param.gNum(2), voxel_param.gNum(1)],sPath+1);
plot([x;x(1)],[y;y(1)],'-.w','MarkerEdgeColor','k', 'LineWidth', 2)
pause(0.1)
    
if ~isempty(walls)
WALLS2{count,1} = walls;
LABELS2{count,1} = label;
PATH2{count,1} = sPath;
sPath(end+1) = sPath(1);
[y,x] = ind2sub([voxel_param.gNum(2), voxel_param.gNum(1)], sPath+1);
mask = mex_fill_inner(x, y, voxel_param.gNum(1), voxel_param.gNum(2));
mask = reshape(mask,voxel_param.gNum(2), voxel_param.gNum(1));
constraintmask(find(mask>0)) = 1;
count = count+1;
end
end



end

