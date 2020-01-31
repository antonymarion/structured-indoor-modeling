%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Room Segmentation
% Input: POINT: point cloud data, CAM: camera index, CAMLIST: camera position
% Ouput: floorplan data
% author: Satoshi Ikehata
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fsc_param = init_freespace_clustering(freespace_evidence, options, voxel_param, numClusters)
% fsc_param = INIT_FREESPACE_CLUSTERING(freespace_evidence, options, voxel_param, numClusters)
% Compute inital medoids and subsampled boundary, region indeces 
% mergin_boundary = options.mergin_boundary/options.leafsize(1);
% mergin_region = options.mergin_region/options.leafsize(1);
% eps_boundary = options.eps_boundary;
% eps_region = options.eps_region;

%% operate boudary freespace
[BW2_resized, BW2] = operate_freespace_evidence(freespace_evidence>0, 1, options.mor_r, 10, options.mw_min_length/options.leafsize(1)); % mor_r = 1 for real data (sparse) 3~5 for synthetic (sparse) datasets

%% extract boundary pixels
binary = BW2_resized;
boundary = mex_findBoundary(voxel_param.gNum(1), voxel_param.gNum(2), double(binary(:)));
boundary_indeces = find(boundary == 1);
boundary_indeces = boundary_indeces(1:end);

%% extract region pixels
[y, x] = meshgrid(1:floor(options.mergin_region/options.leafsize(1)):voxel_param.gNum(2),1:floor(options.mergin_region/options.leafsize(1)):voxel_param.gNum(1));
region_indeces = intersect(find(binary(:) == 1 & boundary(:) == 0), (x(:)-1)*voxel_param.gNum(2)+y(:));

%% subsample boundary & region pixelsr
random_indeces = randperm(length(boundary_indeces));
boundary_indeces = boundary_indeces(random_indeces(1:floor(options.mergin_boundary/options.leafsize(1)):end));
[yBoundary xBoundary] = ind2sub([voxel_param.gNum(2) voxel_param.gNum(1)], boundary_indeces);
[yRegion xRegion] = ind2sub([voxel_param.gNum(2) voxel_param.gNum(1)],region_indeces);

%% compute visibility features
display('Computing visibility features...')
seed = [xRegion yRegion];
target= [xBoundary yBoundary];

feature = mex_FeatureVisibility2D(voxel_param.gNum(1), voxel_param.gNum(2), binary(:), seed, target);
feature = reshape(feature, size(seed,1), size(target,1));

countFeatures = sum(feature);

%% cleaning Boundary and Region
validIndexBoundary = find(countFeatures >= options.eps_boundary);
countFeatures = sum(feature');
validIndexRegion = find(countFeatures >= options.eps_region);

boundary_indeces = boundary_indeces(validIndexBoundary);
region_indeces = region_indeces(validIndexRegion);
[yBoundary, xBoundary] = ind2sub([voxel_param.gNum(2) voxel_param.gNum(1)],boundary_indeces);
[yRegion, xRegion] = ind2sub([voxel_param.gNum(2) voxel_param.gNum(1)], region_indeces);

% add cluster center
display(sprintf('#Region = %d, #Boundary = %d',length(validIndexRegion), length(validIndexBoundary)));
feature = feature(validIndexRegion, validIndexBoundary);

% normalize the visiblity feature
countFeatures = sum(feature);

feature2 = mex_normalize_feature(double(feature), double(countFeatures(:)));
feature2 = reshape(feature2, size(feature,1), size(feature,2));

feature = feature2;

% initialize cluster center
ind = randperm(length(region_indeces));
init_cluster_indeces = ind(1:min(numClusters,floor(length(ind))));

indexRegion = (xRegion-1)*voxel_param.gNum(2) + yRegion;
fsc_param = struct('binary_mask', binary, 'binary_original', BW2, 'init_cluster_indeces', init_cluster_indeces, 'feature', feature, 'xRegion', xRegion, 'yRegion', yRegion, 'indexRegion', indexRegion, 'xBoundary', xBoundary, 'yBoundary', yBoundary);

end

