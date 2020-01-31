%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paper: Structured Indoor Modeling
% published by: IEEE International Conference on Computer Vision (ICCV) 2015
% authors: Satoshi Ikehata, Hang Yan, Yasutaka Furukawa
% contact: sikehata@wustl.edu, yanhang@seas.wustl.edu,
% furukawa@wustl.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function propagated = room_segmentation(freespace_evidence, options, voxel_param)
% propagated = ROOM_ SEGMENTATION(freespace_evidence, options, voxel_param)
% Structure grammar rule: Room segmentation
% Input: Scene Node (represented by 2D free-space evidence)
% Output: Room Nodes (represented by the segmented free-space evidence)

display('[ROOM SEGMENTATION]');
%% parameter setup

numInitClusters = 200; % number of initial clusters
maxIter = 20; % maximum number of kmeans iterations
display(sprintf('#init clusters = %d, max iter = %d', numInitClusters, maxIter));

% initalize cluster centers and samples
fsc_param = init_freespace_clustering(freespace_evidence, options, voxel_param, numInitClusters);

figure('WindowStyle', 'docked', 'Name', 'Region and Boundary'), imagesc(fsc_param.binary_mask), axis equal, hold on
plot(fsc_param.xBoundary, fsc_param.yBoundary, 'or', 'MarkerSize',2)
plot(fsc_param.xRegion, fsc_param.yRegion, '.w', 'MarkerSize',1)
plot(fsc_param.xRegion(fsc_param.init_cluster_indeces),fsc_param.yRegion(fsc_param.init_cluster_indeces), '*g', 'MarkerSize', 10)

%% perform free-space segmentation
display('Running k-medoids clustering...');
[minIndex, minDistance, ambiguity, clusterCenters, indexOnBoundary] = apply_room_segmentation(fsc_param.feature, fsc_param.binary_mask, fsc_param.init_cluster_indeces, fsc_param.indexRegion, maxIter, options.pAssign, options.pMerge, fsc_param.xBoundary, fsc_param.yBoundary);

minIndexMap = zeros(voxel_param.gNum(1)*voxel_param.gNum(2), 1);
minDistanceMap = zeros(voxel_param.gNum(1)*voxel_param.gNum(2), 1);
ambiguityMap = zeros(voxel_param.gNum(1)*voxel_param.gNum(2), 1);

minIndexMap(fsc_param.indexRegion) = minIndex;
minIndexMap = reshape(minIndexMap, voxel_param.gNum(2),voxel_param.gNum(1));
minDistanceMap(fsc_param.indexRegion) = minDistance;
minDistanceMap = reshape(minDistanceMap, voxel_param.gNum(2), voxel_param.gNum(1));
ambiguityMap(fsc_param.indexRegion) = ambiguity;
ambiguityMap = reshape(ambiguityMap, voxel_param.gNum(2), voxel_param.gNum(1));


%% output index and distance maps
minIndexMap_original = zeros(voxel_param.gNum(2),voxel_param.gNum(1));
minDistanceMap_original = zeros(voxel_param.gNum(2),voxel_param.gNum(1));
ambiguityMap_original = zeros(voxel_param.gNum(2),voxel_param.gNum(1));

[X Y] = meshgrid([1:voxel_param.gNum(1)], [1:voxel_param.gNum(2)]);
index = (X(:)-1)*voxel_param.gNum(2) + Y(:);

minIndexMap_original(index) = minIndexMap(:);
minIndexMap_original = reshape(minIndexMap_original, voxel_param.gNum(2),voxel_param.gNum(1));
minDistanceMap_original(index) = minDistanceMap(:);
minDistanceMap_original = reshape(minDistanceMap_original, voxel_param.gNum(2),voxel_param.gNum(1));
ambiguityMap_original(index) = ambiguityMap(:);
ambiguityMap_original = reshape(ambiguityMap_original, voxel_param.gNum(2),voxel_param.gNum(1));

[cluster_fs, distance_fs, ambiguity_fs, M] = upsample_index_distance_maps(fsc_param.binary_original, minIndexMap_original, minDistanceMap_original, ambiguityMap_original);
propagated = refine_indexmap(M, fsc_param.binary_original, cluster_fs);

display(sprintf('Result of room segmentation rule: 1 scene -> %d room nodes', length(unique(propagated))-1));
end