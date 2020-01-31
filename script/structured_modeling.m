%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paper: Structured Indoor Modeling
% published by: IEEE International Conference on Computer Vision (ICCV) 2015
% authors: Satoshi Ikehata, Hang Yan, Yasutaka Furukawa
% contact: sikehata@wustl.edu, yanhang@seas.wustl.edu,
% furukawa@wustl.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% In this implementation, we sequentailly apply structure grammer rules as
% follow
% 0. Preprocess
% 1. Room segmentation
% 2. Room reconstruction (wall)
% 3. Room addition + reconstruction (wall)
% 4. Room merging and door addition (1st)
% 5. Room merging and door addition (2nd) 
% 6. Room reconstruction (floor, ceiling)
% 7. Export structure graph (only contains basic architectural elements)
% 8. Detail reconstruction



%% (0-a) Noisy point removal and Manhattan-World Coordinate Estimation
[options, core_point_idx] = point_noise_removal(POINT, options); % if the input point cloud is noisy perform this noise removal (you can replace here by any other algorithms to get clean point cloud)
[POINT, CAMLIST, manhattanworld_rotation] = poc_manhattan(POINT, CAMLIST, core_point_idx);


%% (0-b) Compute the evidences (2D/3D point/free-space evidences) 
[Evidences, voxel_param] = compute_evidences(POINT, CAM, CAMLIST, core_point_idx, options);

figure('WindowStyle', 'docked', 'Name', 'freespace-evidence (dense)'); imagesc(Evidences.freespace_evidence), axis equal % save intermediate images
saveas(gcf,sprintf('%s/freespace2d.png',intermediate_images), 'png');
figure('WindowStyle', 'docked', 'Name', 'point-evidence (all)'); imagesc(Evidences.point_evidence), axis equal
saveas(gcf,sprintf('%s/point2d.png',intermediate_images), 'png');
figure('WindowStyle', 'docked', 'Name', 'point-evidence (all)'); imagesc(Evidences.point_evidence_wall), axis equal
saveas(gcf,sprintf('%s/wall2d.png',intermediate_images), 'png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Run structured modeling %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% (1) Apply room segmentation rule (a scene node -> room nodes)
propagated = room_segmentation(Evidences.freespace_evidence, options, voxel_param);
Evidences.binarymask = double(propagated>0); 
figure('WindowStyle', 'docked', 'Name', 'K-medois segmentation'), imagesc(propagated), axis equal
saveas(gcf,sprintf('%s/kmedoids.png',intermediate_images), 'png');


%% (2) Apply room reconstruction rule (room nodes->walls)
[labels, connectivity]  = find_connectivity(propagated);
[WALLS, LABELS, PATH] = room_reconstruction_wall(setdiff(unique(propagated),0),  zeros(voxel_param.gNum(2), voxel_param.gNum(1)), propagated, Evidences, connectivity, wall_param, voxel_param, options);
filename = sprintf('%s/graph_roomreconstruction.png',intermediate_images);
graphshow(voxel_param, WALLS, filename, []);


%% (3) Appy room addition rule (with room reconstruction for new room nodes)
[WALLS_WITH_UNEXPLAINED, LABELS_WITH_UNEXPLAINED, PATH_WITH_UNEXPLAINED, unexplained_indeces] = room_addition_reconstruction(WALLS, LABELS, PATH, Evidences, wall_param, voxel_param, options);
filename = sprintf('%s/graph_roomaddition.png',intermediate_images);
graphshow(voxel_param, WALLS_WITH_UNEXPLAINED, filename, []);


%% (4) Apply door addition and room merging rules (1st)
[mergeList, doorList] = room_connection_analysis(WALLS_WITH_UNEXPLAINED,  LABELS_WITH_UNEXPLAINED, PATH_WITH_UNEXPLAINED, Evidences, voxel_param, door_param, options); % Check graph condition to find if either of rules could be applied (1st)
if isempty(doorList)
    return;
end

[WALLS_WITH_UNEXPLAINED_REMOVED, LABELS_WITH_UNEXPLAINED_REMOVED, PATH_WITH_UNEXPLAINED_REMOVED, doorList_removed, mergeList_removed] = delete_empty_rooms(WALLS_WITH_UNEXPLAINED, LABELS_WITH_UNEXPLAINED, PATH_WITH_UNEXPLAINED, voxel_param, doorList, unexplained_indeces); % Remove room nodes with no connection to other room nodes
filename = sprintf('%s/connection_analysis_1st.png',intermediate_images);
graphshow(voxel_param, WALLS_WITH_UNEXPLAINED_REMOVED, filename, doorList);


if ~isempty(mergeList)
[WALLS_MERGED, LABELS_MERGED, PATH_MERGED] = merge_rooms(mergeList_removed, doorList_removed, WALLS_WITH_UNEXPLAINED_REMOVED,  LABELS_WITH_UNEXPLAINED_REMOVED, PATH_WITH_UNEXPLAINED_REMOVED, Evidences, wall_param, voxel_param, options);
[WALLS_MERGED, LABELS_MERGED, PATH_MERGED] = delete_room_in_room(WALLS_MERGED, LABELS_MERGED, PATH_MERGED, voxel_param);

%% (5) Apply door addition and room merging rules (2nd)
[mergeList, doorList] = room_connection_analysis(WALLS_MERGED, LABELS_MERGED, PATH_MERGED, Evidences, voxel_param, door_param, options); % Check graph condition to find if either of rules could be applied (2nd)
filename = sprintf('%s/connection_analysis_2nd.png',intermediate_images);
graphshow(voxel_param, WALLS_MERGED, filename, doorList);
else
WALLS_MERGED = WALLS_WITH_UNEXPLAINED;
LABELS_MERGED = LABELS_WITH_UNEXPLAINED;
PATH_MERGED = PATH_WITH_UNEXPLAINED;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this implementation, we stop the application of room merging and door
% addition after ONE iteration. It means all MERGE connection at the 2nd connection analysis are discarded. 
%You can change this part to apply these grammaer rules untile converges (may not improve the result, however)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1: length(doorList), doorList{i}.aisleFlag = 0; end
[WALLS_REMOVED, LABELS_REMOVED, PATH_REMOVED, doorList_removed] = delete_empty_rooms(WALLS_MERGED, LABELS_MERGED, PATH_MERGED, voxel_param, doorList, []);
propagated2 = geodesic_propagation(path_to_mask(PATH_REMOVED, voxel_param), double(Evidences.freespace_evidence>0));
filename = sprintf('%s/connection_analysis_final.png',intermediate_images);
graphshow(voxel_param, WALLS_REMOVED, filename, doorList_removed);


%% (6) Apply room reconstruction rule (floor, ceiling)
[POINT_room, CAM_room] = assign_room_point(POINT, CAM, CAMLIST, options, voxel_param, core_point_idx, propagated2);
[floorHeight, ceilHeight] = floor_ceil_reconstruction(POINT_room, CAM_room, CAMLIST);
floorHeight(:) = min(floorHeight); % assign the same floor height to all room nodes (remove this allows different floor heights)

display('Completed the structure graph recovery of basic architectural elements');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Before here, the structure graph about the basic architectural elements are recovered 
% Details and objects are still missing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% (7) Export structure graph (only basic architectural elements)
doorListFinal = cell(1,1);
count = 0;
for i = 1 : length(doorList_removed)
    if doorList_removed{i}.aisleFlag == 0
        doorListFinal{count+1} = doorList_removed{i};
        count = count + 1;
    end
end
filename = sprintf('%s/floorplan.txt', geometry_dir);
write_structure_graph_basic(filename, WALLS_REMOVED, LABELS_REMOVED, doorListFinal, floorHeight, ceilHeight, manhattanworld_rotation);
convert_floorplan_to_ply(filename);

%% 


