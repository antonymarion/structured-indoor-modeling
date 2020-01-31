%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paper: Structured Indoor Modeling
% published by: IEEE International Conference on Computer Vision (ICCV) 2015
% authors: Satoshi Ikehata, Hang Yan, Yasutaka Furukawa
% contact: sikehata@wustl.edu, yanhang@seas.wustl.edu,
% furukawa@wustl.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [WALLS2, LABELS2, PATH2, unexplained_index] = room_addition_reconstruction(WALLS, LABELS, PATH, Evidences, wall_param, voxel_param, options)
% [WALLS2, LABELS2, PATH2, unexplained_index] = ROOM_ADDITION_RECONSTRUCTION(WALLS, LABELS, PATH, Evidences, wall_param, voxel_param, options)
% Apply room addition rule to add a new room node and reconstruct walls via
% room reconstruction rule 

display('[ROOM ADDITION & RECONSTRUCTION]');
    display('Searching for unexplained freespace regions...');
BW = double(Evidences.binarymask);
newClusterMap = path_to_mask(PATH, voxel_param);
propagated = geodesic_propagation(newClusterMap, BW);

indeces = setdiff(unique(newClusterMap),0);
area = zeros(length(indeces), 1);
for i = 1 : length(indeces)
   area(i) = length(find(newClusterMap==indeces(i)));
end
thresh1 = 0.05*mean(area(:));% default 0.05
thresh2 = 0.5; % default 0.5 

WALLS2 = WALLS;
LABELS2 = LABELS;
PATH2 = PATH;


% process on unexplained region
[walls, labels, path] = reconstruct_additional_room(Evidences, newClusterMap, propagated, thresh1, thresh2, wall_param, voxel_param, options);
unexplained_index = [];
while(1)
if ~isempty(cell2mat(walls))
    unexplained_index = [unexplained_index;[(1:length(walls))+length(WALLS2)]'];
    WALLS2 = [WALLS2;walls]; LABELS2 = [LABELS2;labels]; PATH2 = [PATH2;path];
    [WALLS2, LABELS2, PATH2, newClusterMap] = arrange_room_nodes(WALLS2, LABELS2, PATH2, voxel_param);
    propagated = geodesic_propagation(newClusterMap, BW);
    display('Searching for unexplained freespace regions...');
    [walls, labels, path] = reconstruct_additional_room(Evidences, newClusterMap, propagated, thresh1, thresh2, wall_param, voxel_param, options);
else  
    break;
end
end

display(sprintf('Updated structure graph. # Room %d -> %d', length(WALLS), length(WALLS2)));

end
