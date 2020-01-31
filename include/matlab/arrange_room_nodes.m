%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paper: Structured Indoor Modeling
% published by: IEEE International Conference on Computer Vision (ICCV) 2015
% authors: Satoshi Ikehata, Hang Yan, Yasutaka Furukawa
% contact: sikehata@wustl.edu, yanhang@seas.wustl.edu,
% furukawa@wustl.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [WALLS, LABELS, PATH, newClusterMap] = arrange_room_nodes(WALLS, LABELS, PATH, voxel_param)
% [WALLS, LABELS, PATH, newClusterMap] = ARRANGE_ROOM_NODES(WALLS, LABELS, PATH, voxel_param)
% Update all nodes after the graph update

valid_cluster_list = [];
count = 0;

newClusterMap = zeros(voxel_param.gNum(2), voxel_param.gNum(1));
count = 0;
for i = 1 : length(PATH)
    if ~isempty(PATH{i})
    count = count + 1;
    valid_cluster_list(count) = i;
    sPath = PATH{i};
    sPath(end+1) = sPath(1);
    labels = LABELS{i};
    [y x] = ind2sub([voxel_param.gNum(2), voxel_param.gNum(1)],sPath+1); % should be the same size
    mask = mex_fill_inner(x, y, voxel_param.gNum(1), voxel_param.gNum(2));
    mask = reshape(mask,voxel_param.gNum(2), voxel_param.gNum(1));
    newClusterMap(find(mask==1)) = i;
    end
end

newClusterMap_ = zeros(voxel_param.gNum(2), voxel_param.gNum(1));
WALLS_ = cell(length(valid_cluster_list),1);
LABELS_ = cell(length(valid_cluster_list),1);
PATH_ = cell(length(valid_cluster_list),1);

for i = 1 : length(valid_cluster_list)
    newClusterMap_(find(newClusterMap == valid_cluster_list(i))) = i;   
    WALLS_{i} = WALLS{valid_cluster_list(i)};
    LABELS_{i} = LABELS{valid_cluster_list(i)};
    PATH_{i} = PATH{valid_cluster_list(i)};
end
newClusterMap = newClusterMap_;
WALLS = WALLS_;
LABELS = LABELS_;
PATH = PATH_;
end

