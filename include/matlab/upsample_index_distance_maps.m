%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paper: Structured Indoor Modeling
% published by: IEEE International Conference on Computer Vision (ICCV) 2015
% authors: Satoshi Ikehata, Hang Yan, Yasutaka Furukawa
% contact: sikehata@wustl.edu, yanhang@seas.wustl.edu,
% furukawa@wustl.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cluster_fs, distance_fs, ambiguity_fs, M] = upsample_index_distance_maps(BW, minIndexMap_original, minDistanceMap_original, ambiguityMap_original)
% [cluster_fs, distance_fs, ambiguity_fs, M] = UPSAMPLE_INDEX_DISTANCE_MAPS(BW, minIndexMap_original, minDistanceMap_original, ambiguityMap_original)
% Upsample the subsampled segmented free-space evidence upon the original
% resolution using geodesic upsamling

%% propagate to free-space evidence using geodesic filtering transform
[h,w] = size(BW);
BW_copy = double(BW);
BW_copy(find(BW_copy==0)) = rand(size(find(BW_copy==0)));
[dIf dIb idf idb] = calc_dIandId_gray(BW_copy(:), h, w);

%% for Index Map
idzero = find(minIndexMap_original>0);
L0 = [1:h*w]';
L0 = reshape(L0,h,w);
[X, Y] = meshgrid([1:w],[1:h]);
id = (X(:)-1)*h + Y(:);
M0 = 1.0e12*ones(h*w,1);
M0(idzero) = 0;
M0 = reshape(M0,h,w);

% geodesic distance transform
r = 1;
lambda = 1.0e12;
[M, L] = mex_gdt(M0, L0, dIf, dIb, r, lambda, 2);
M = reshape(M,h,w);
L = reshape(L,h,w);

cluster_fs = minIndexMap_original(L);
cluster_fs(BW==0) = 0;


%% for Distance Map
idzero = find(minDistanceMap_original>0);
L0 = [1:h*w]';
L0 = reshape(L0,h,w);
[X, Y] = meshgrid([1:w],[1:h]);
M0 = 1.0e12*ones(h*w,1);
M0(idzero) = 0;
M0 = reshape(M0,h,w);

% geodesic distance transform
r = 1;
lambda = 1.0e12;
[M, L] = mex_gdt(M0, L0, dIf, dIb, r, lambda, 2);
M = reshape(M,h,w);
L = reshape(L,h,w);

distance_fs = minDistanceMap_original(L);
distance_fs(BW==0) = 0;

%% for Ambiguity Map
idzero = find(ambiguityMap_original>0);
L0 = [1:h*w]';
L0 = reshape(L0,h,w);
[X, Y] = meshgrid([1:w],[1:h]);
M0 = 1.0e12*ones(h*w,1);
M0(idzero) = 0;
M0 = reshape(M0,h,w);
% geodesic distance transform

r = 1;
lambda = 1.0e12;
[M, L] = mex_gdt(M0, L0, dIf, dIb, r, lambda, 2);
M = reshape(M,h,w);
L = reshape(L,h,w);

ambiguity_fs = ambiguityMap_original(L);
ambiguity_fs(BW==0) = 0;
end

