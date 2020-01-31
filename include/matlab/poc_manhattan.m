function [POINT, CAMLIST, manhattanworld_rotation] = poc_manhattan(POINT, CAMLIST, core_point_idx)
% [POINT, CAMLIST, manhattanworld_rotation] = POC_MANHATTAN(POINT, CAMLIST, core_point_idx)
% Compute Manhattan-World coordinate axes and rotate points along their
% directions

display('...computing Manhattan-World coordinate')
point = cell2mat(POINT);
random_indeces = randperm(length(core_point_idx)); 
% random_indeces = random_indeces(1:floor(length(core_point_idx)/100000):end);
random_indeces = core_point_idx(random_indeces(1:10:end)); % only use every 10 points for the efficiency
manhattanworld_rotation = find_manhattan(point(random_indeces,3:5), point(random_indeces,9:11));

% visualization optional
pocshow((manhattanworld_rotation*point(random_indeces,3:5)')', -1);
%% rotate points
CAMLIST = (manhattanworld_rotation*CAMLIST')';
POINT = rotate_point_3d(POINT, manhattanworld_rotation);

end

