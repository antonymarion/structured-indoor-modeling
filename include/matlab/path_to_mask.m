%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paper: Structured Indoor Modeling
% published by: IEEE International Conference on Computer Vision (ICCV) 2015
% authors: Satoshi Ikehata, Hang Yan, Yasutaka Furukawa
% contact: sikehata@wustl.edu, yanhang@seas.wustl.edu,
% furukawa@wustl.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function newClusterMap = path_to_mask(PATH, voxel_param)
% newClusterMap = PATH_TO_MASK(PATH, voxel_param)
% Convert 1D wall path to a polygon in the 2-D domain
newClusterMap = zeros(voxel_param.gNum(2), voxel_param.gNum(1));
for i = 1 : length(PATH)
    if ~isempty(PATH{i})
    sPath = PATH{i};
    sPath(end+1) = sPath(1);
    [y x] = ind2sub([voxel_param.gNum(2), voxel_param.gNum(1)],sPath+1); % should be the same size
    mask = mex_fill_inner(x, y, voxel_param.gNum(1), voxel_param.gNum(2));
    mask = reshape(mask,voxel_param.gNum(2), voxel_param.gNum(1));
    newClusterMap(find(mask==1)) = i;
    end
end
end

