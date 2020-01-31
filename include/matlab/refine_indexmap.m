%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paper: Structured Indoor Modeling
% published by: IEEE International Conference on Computer Vision (ICCV) 2015
% authors: Satoshi Ikehata, Hang Yan, Yasutaka Furukawa
% contact: sikehata@wustl.edu, yanhang@seas.wustl.edu,
% furukawa@wustl.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cluster_fs_refine2 = refine_indexmap(M, BW, cluster_fs)
% cluster_fs_refine2 = REFINE_INDEXMAP(M, BW, cluster_fs)
% remove isolated cluster (each room node should contain connected region)

[h,w] = size(BW);
%% propagate to free-space evidence using geodesic filtering transform
cluster_fs(M>0.9e10) = 0;

cluster_fs_refine = delete_isolated_cluster(cluster_fs);

BW_copy = double(BW);
BW_copy(find(BW_copy==0)) = rand(size(find(BW_copy==0)));
[dIf dIb idf idb] = calc_dIandId_gray(BW_copy(:), h, w);

%% for Index Map
M = zeros(h*w,1);
L = zeros(h*w,1);
idzero = find(cluster_fs_refine>0);
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
cluster_fs_refine2 = cluster_fs_refine(L);
cluster_fs_refine2(BW==0) = 0;
end

