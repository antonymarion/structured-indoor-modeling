%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paper: Structured Indoor Modeling
% published by: IEEE International Conference on Computer Vision (ICCV) 2015
% authors: Satoshi Ikehata, Hang Yan, Yasutaka Furukawa
% contact: sikehata@wustl.edu, yanhang@seas.wustl.edu,
% furukawa@wustl.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [WALLS, LABELS, PATH] = room_reconstruction_wall(core_indeces, mask, cluster_map, Evidences, connectivity, wall_param, voxel_param, options)
% [WALLS, LABELS, PATH] = ROOM_RECONSTRUCTION_WALL(core_indeces, mask, cluster_map, Evidences, connectivity, wall_param, voxel_param, options)
% Recover the 1D wall path for each room reagion

display('[ROOM RECONSTRUCTION (WALL)]');
display(sprintf('# of room nodes = %d', length(connectivity)));
numcluster = length(core_indeces);
% core_indeces = setdiff(unique(cluster_map),[target_index(:);0]);
num_connectivity = zeros(numcluster, 1);
for i = 1 : numcluster
    num_connectivity(i) = length(connectivity{core_indeces(i)});
end

order_clusters = sortrows([[1:numcluster]',num_connectivity], 2);
constraintmask = mask;

%% computing the core freespace evidence
display('Computing the core free-space evidence...')

coremask = zeros(voxel_param.gNum(2),voxel_param.gNum(1));
coremask_large  = zeros(voxel_param.gNum(2),voxel_param.gNum(1));
corecorners = zeros(length(core_indeces),4);

for i = 1 : length(core_indeces)
    index_cluster = core_indeces(order_clusters(i,1));
    supportmask = cluster_map == index_cluster;
    supportmask(constraintmask>0) = 0;
    inv_mask = 1-supportmask;    
    [yy xx] = ind2sub([voxel_param.gNum(2), voxel_param.gNum(1)], find(supportmask>0));
    integral_inv_mask = integralImage(inv_mask);    
    integral_mask_pt = integralImage(Evidences.point_evidence);
    X = mex_compute_core_freespace(double(supportmask), double(integral_inv_mask), double(integral_mask_pt), min(xx)-1, max(xx)-1, min(yy)-1, max(yy)-1);
    corecorners(i, :) = X;
    [XX YY] = meshgrid([X(1)+wall_param.core_freespace_margin:X(3)- wall_param.core_freespace_margin],  [X(2)+wall_param.core_freespace_margin:X(4)- wall_param.core_freespace_margin]);    
    [XX2 YY2] = meshgrid([X(1):X(3)],  [X(2):X(4)]);    
    index = (XX-1)*voxel_param.gNum(2) + YY;
    index_large = (XX2-1)*voxel_param.gNum(2) + YY2;
    coremask(index) = index_cluster;   
    coremask_large(index_large) = index_cluster;
    display(sprintf('...Size of CFE (room %d) = %d', i, length(index)));
end

%% computing the shortest path
WALLS = cell(1,1);
LABELS = cell(1,1);
PATH = cell(1,1);

count = 0;
for i = 1 : length(core_indeces)
    display(sprintf('Room %d: Wall nodes reconstruction (# room node is %d)', i, length(core_indeces)));
    index_cluster = core_indeces(order_clusters(i,1));
    invalidmask = double(coremask_large ~= index_cluster & coremask_large > 0);    
    cmask = cluster_map > 0;
    cmask = imdilate(cmask,strel('disk',20));
    constraintmask(cmask==0) = 1;
    point_evidence_cluster = zeros(size(Evidences.point_evidence_wall));
    point_evidence_cluster(cluster_map==index_cluster) = Evidences.point_evidence_wall(cluster_map == index_cluster);
    

    
    %% find start-end-line, and set core mask
    
    core_room = double(coremask==index_cluster);
    if isempty(find(core_room==1))
       continue; 
    end
    
    % convert core corners to the points 
    points = cell(4,1);
    corners = zeros(5,2);
    corners(1,:) = [corecorners(i,1), corecorners(i,2)];
    corners(2,:) = [corecorners(i,3), corecorners(i,2)];
    corners(3,:) = [corecorners(i,3), corecorners(i,4)];
    corners(4,:) = [corecorners(i,1), corecorners(i,4)];
    corners(5,:) = corners(1,:);
    
    for t = 1 : size(corners,1)-1
    d = [corners(t+1,1)-corners(t,1),corners(t+1,2)-corners(t,2)];
    numpix = norm(d);
    d = d/numpix;
    points{t,1} = [corners(t,1)+([1:numpix+1]-1)*d(1);corners(t,2)+([1:numpix+1]-1)*d(2)]';
    end      
    points = cell2mat(points);   
    path_evidence = zeros(size(points,1),1);
    [xStart, xGoal, core_room] = mex_find_start_end(core_room, points'-1,  double(path_evidence), double(point_evidence_cluster), double(zeros(size(point_evidence_cluster))));
    core_room(constraintmask==1) = 1;
    core_room(invalidmask==1) = 1;    
    display('... Running 1st Dijkstra algorithm')
    
%     figure, imagesc(reshape(core_room, voxel_param.gNum(2), voxel_param.gNum(1))), axis equal, hold on
%     plot(xStart(1)+1, xStart(2)+1, 'wo')
%     plot(xGoal(1)+1, xGoal(2)+1, 'wo')
    %% path reconstruction
%     tic;
    [walls, label, sPath] = room_wall_path_extraction_normal(xStart, xGoal, point_evidence_cluster, Evidences.freespace_evidence_sparse, wall_param, voxel_param, options, core_room, Evidences.nx, Evidences.ny);   
%     time = toc;
%     display(sprintf('elapsed time is %.02f sec', time));
        
    %% path evaluation and find new start-end point
    [y x] = ind2sub([voxel_param.gNum(2), voxel_param.gNum(1)], [sPath(:)+1;sPath(1)+1]);    
    points = cell(length(sPath), 1);
    path_evidence = cell(length(sPath), 1);

    for t = 1 : length(sPath)
        d = [x(t+1)-x(t),y(t+1)-y(t)];
        numpix = norm(d)+1;
        d = d/(numpix-1);
        X = [x(t)+([1:numpix]-1)*d(1);y(t)+([1:numpix]-1)*d(2)]';
        index_path = (X(:,1)-1)*voxel_param.gNum(2) + X(:,2);
        [yy xx] = ind2sub([voxel_param.gNum(2), voxel_param.gNum(1)], index_path);
        points{t,1} = [xx,yy];
        path_evidence{t,1} = sum(Evidences.point_evidence_wall(index_path))*ones(length(xx),1);
    end
 
   
   
    points = cell2mat(points);
    path_evidence = cell2mat(path_evidence);
    core_room = double(coremask==index_cluster);    
    
    % The first start-end line may not be optimal. 
    % Update it after the first path reconstruction.
    [xStart, xGoal, core_room] = mex_find_start_end(core_room, points'-1,  double(path_evidence), double(point_evidence_cluster),double(zeros(size(point_evidence_cluster)))); 
    core_room(constraintmask==1) = 1;
    core_room(invalidmask==1) = 1;    
    
    display('... Running 2nd Dijkstra algorithm')

    [walls, label, sPath] = room_wall_path_extraction_normal(xStart, xGoal, point_evidence_cluster, Evidences.freespace_evidence_sparse, wall_param, voxel_param, options, core_room, Evidences.nx, Evidences.ny);   

           
    % perform filtering 
    kernel = fspecial('gaussian', 2*ceil(wall_param.margin_size/options.leafsize(1))+1, wall_param.margin_sigma);
    test = Evidences.freespace_evidence_sparse;
    figure('WindowStyle', 'docked', 'Name', sprintf('Room %d',i)), imagesc(reshape(test, voxel_param.gNum(2), voxel_param.gNum(1))), axis equal
    hold on
    plot(xStart(1)+1,xStart(2)+1,'or','MarkerEdgeColor','g','MarkerSize',6)
    plot(xGoal(1)+1,xGoal(2)+1,'*r','MarkerEdgeColor','b','MarkerSize',6)
    [y x] = ind2sub([voxel_param.gNum(2), voxel_param.gNum(1)],sPath+1);
    plot([x;x(1)],[y;y(1)],'-.w','MarkerEdgeColor','k', 'LineWidth', 2)
    pause(0.1)
    
    
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
    count = count+1;
    WALLS{count,1} = walls;
    LABELS{count,1} = label;
    PATH{count,1} = sPath;
    sPath(end+1) = sPath(1);
    [y,x] = ind2sub([voxel_param.gNum(2), voxel_param.gNum(1)], sPath+1);
    mask = mex_fill_inner(x, y, voxel_param.gNum(1), voxel_param.gNum(2));
    mask = reshape(mask,voxel_param.gNum(2), voxel_param.gNum(1));    
    constraintmask(find(mask>0)) = 1;
    display(sprintf('# of wall nodes is %d', length(sPath)));
    end  
end

end