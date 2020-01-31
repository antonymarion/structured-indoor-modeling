%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paper: Structured Indoor Modeling
% published by: IEEE International Conference on Computer Vision (ICCV) 2015
% authors: Satoshi Ikehata, Hang Yan, Yasutaka Furukawa
% contact: sikehata@wustl.edu, yanhang@seas.wustl.edu,
% furukawa@wustl.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [WALLS, LABELS, PATH] = reconstruct_additional_room(Evidences, core_cluster_map, propagated, thresh1, thresh2, wall_param, voxel_param, options)
% [WALLS, LABELS, PATH] = RECONSTRUCT_ADDITIONAL_ROOM(Evidences, core_cluster_map, propagated, thresh1, thresh2, wall_param, voxel_param, options)
% Apply room reconstruction rule for an additional rooms
  

    [h, w] = size(core_cluster_map);
    assigned_map = core_cluster_map == 0;
    unassigned_map = assigned_map & propagated > 0;
    CC = bwconncomp(unassigned_map);
    count = 0;
    
    WALLS = cell(1,1);
    LABELS = cell(1,1);
    PATH = cell(1,1);

    constraintmask = zeros(voxel_param.gNum(2), voxel_param.gNum(1));
    new_core_cluster_map = core_cluster_map;
    
    display(sprintf('Found %d unexplained regions...', CC.NumObjects));
    for i = 1 : CC.NumObjects

       region_indeces = CC.PixelIdxList{i};
       label_list = unique(propagated(region_indeces));
       label_list = setdiff(label_list, 0);
       fs = Evidences.freespace_evidence > 0;
       pt = Evidences.point_evidence > 0;
       
       if(length(label_list) > 1 & length(region_indeces) > thresh1 & sum(pt(region_indeces))/sum(fs(region_indeces)) > thresh2)     
           
    
     
        [h, w] = size(propagated);
        mask = zeros(h, w);
        mask(region_indeces) = propagated(region_indeces);        
        se = strel('disk',3);
%         figure, imagesc(mask), axis equal
        temp = imerode(mask>0, se);
        temp = imdilate(temp>0, se);
%      figure, imagesc(mask), axis equal
        mask(temp==0) = 0;
        supportmask = mask>0;
        inv_mask = 1-supportmask;    
        [yy xx] = ind2sub([h, w], find(mask>0));
        integral_inv_mask = integralImage(inv_mask);    
%         tic;
        X = mex_compute_core_freespace_unexplained(double(supportmask), double(integral_inv_mask), mask, min(xx)-1, max(xx)-1, min(yy)-1, max(yy)-1); % X is in matlab coordinate (+1)
%         toc;
        corecorners(i, :) = X;
        [XX YY] = meshgrid([X(1)+wall_param.core_freespace_margin:X(3)- wall_param.core_freespace_margin],  [X(2)+wall_param.core_freespace_margin:X(4)- wall_param.core_freespace_margin]);    
        [XX2 YY2] = meshgrid([X(1):X(3)],  [X(2):X(4)]);    
        index = (XX-1)*voxel_param.gNum(2) + YY;
        index_large = (XX2-1)*voxel_param.gNum(2) + YY2;
        coremask = zeros(voxel_param.gNum(2), voxel_param.gNum(1));
        coremask_large = zeros(voxel_param.gNum(2), voxel_param.gNum(1));
        coremask(index) = 1;   
        coremask_large(index_large) = 1;     
        

        if isempty(index) | length(index_large(:)) < thresh1
            continue;
        end
        
        newIndex = max(unique(new_core_cluster_map))+1;
        new_core_cluster_map(coremask == 1) = newIndex;
        
        invalidmask = core_cluster_map > 0;    
        cmask = propagated > 0;
        cmask = imdilate(cmask,strel('disk',20));
        constraintmask(cmask==0) = 1;
        point_evidence_cluster = zeros(size(Evidences.point_evidence_wall));
        
        mask = mask > 0;
        mask = imopen(mask, strel('disk', 5));
        point_evidence_cluster(mask>0) = Evidences.point_evidence_wall(mask>0);
%         point_evidence_cluster(new_propagated==newIndex) = Evidences.point_evidence_wall(new_propagated==newIndex); % which one is better?
    
        % convert core corners to the points 
        points = cell(4,1);
        corners = zeros(5,2);
        corners(1,:) = [corecorners(i,1), corecorners(i,2)];
        corners(2,:) = [corecorners(i,3), corecorners(i,2)];
        corners(3,:) = [corecorners(i,3), corecorners(i,4)];
        corners(4,:) = [corecorners(i,1), corecorners(i,4)];
        corners(5,:) = corners(1,:);
    
        for i = 1 : size(corners,1)-1
        d = [corners(i+1,1)-corners(i,1),corners(i+1,2)-corners(i,2)];
        numpix = norm(d);
        d = d/numpix;
        points{i,1} = [corners(i,1)+([1:numpix+1]-1)*d(1);corners(i,2)+([1:numpix+1]-1)*d(2)]';
        end      
        points = cell2mat(points);   
        path_evidence = zeros(size(points,1),1);
        [xStart, xGoal, core_room] = mex_find_start_end(coremask, points'-1,  double(path_evidence), double(point_evidence_cluster),double(zeros(size(point_evidence_cluster))));
        core_room(constraintmask==1) = 1;
        core_room(invalidmask==1) = 1;    
        
        display(sprintf('Size of largest unexplained region = %d : Wall nodes reconstruction', abs((X(1)-X(3))*(X(2)-X(4)))));
        display('... Running 1st Dijkstra algorithm')
    
        %% path reconstruction
%         figure, imagesc(reshape(core_room, voxel_param.gNum(2), voxel_param.gNum(1))), axis equal, hold on
%         plot(xStart(1)+1,xStart(2)+1, 'ow')
%         plot(xGoal(1)+1,xGoal(2)+1, '*w')
        
        tic;
        [walls, label, sPath] = room_wall_path_extraction_normal(xStart, xGoal, point_evidence_cluster, Evidences.freespace_evidence_sparse, wall_param, voxel_param, options, core_room, Evidences.nx, Evidences.ny);   
%         time = toc;
%         display(sprintf('elapsed time is %.02f sec', time));
    
        %% path evaluation and find new start-end point
        [y x] = ind2sub([voxel_param.gNum(2), voxel_param.gNum(1)], [sPath(:)+1;sPath(1)+1]);    
        points = cell(length(sPath), 1);
        path_evidence = cell(length(sPath), 1);

        for i = 1 : length(sPath)
            d = [x(i+1)-x(i),y(i+1)-y(i)];
            numpix = norm(d);
            d = d/numpix;
            X = [x(i)+([1:numpix+1]-1)*d(1);y(i)+([1:numpix+1]-1)*d(2)]';
            index_path = (X(:,1)-1)*voxel_param.gNum(2) + X(:,2);
            [yy xx] = ind2sub([voxel_param.gNum(2), voxel_param.gNum(1)], index_path);
            points{i,1} = [xx,yy];
            path_evidence{i,1} = sum(Evidences.point_evidence_wall(index_path))*ones(length(xx),1);
        end

        points = cell2mat(points);
        path_evidence = cell2mat(path_evidence); 

        [xStart, xGoal, core_room] = mex_find_start_end(coremask, points'-1,  double(path_evidence), double(point_evidence_cluster),double(zeros(size(point_evidence_cluster))));
        core_room(constraintmask==1) = 1;
        core_room(invalidmask==1) = 1;    

        display('... Running 2nd Dijkstra algorithm')
%         tic;
        [walls, label, sPath] = room_wall_path_extraction_normal(xStart, xGoal, point_evidence_cluster, Evidences.freespace_evidence_sparse, wall_param, voxel_param, options, core_room, Evidences.nx, Evidences.ny);   
%         time = toc;
%         display(sprintf('elapsed time is %.02f sec', time));  

        test =point_evidence_cluster;
        test(find(core_room==1)) = -2;
        test = test + Evidences.freespace_evidence;
%         figure('WindowStyle', 'docked', 'Name', sprintf('Unexplained: %d', count)), imagesc(reshape(test, voxel_param.gNum(2), voxel_param.gNum(1))), axis equal
%         hold on
%         plot(xStart(1)+1,xStart(2)+1,'or','MarkerEdgeColor','g','MarkerSize',6)
%         plot(xGoal(1)+1,xGoal(2)+1,'*r','MarkerEdgeColor','b','MarkerSize',6)
%         [y x] = ind2sub([voxel_param.gNum(2), voxel_param.gNum(1)],sPath+1);
%         plot([x;x(1)],[y;y(1)],'-.w','MarkerEdgeColor','k', 'LineWidth', 2)
%         pause(0.1)

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
        end            
           
          
       end
    end
end

