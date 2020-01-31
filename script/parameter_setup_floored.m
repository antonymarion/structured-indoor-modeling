% Parameters setup for our floored dataset (excluding office2)

% General parameters
% leafsize: size of a boxel grid in mm
% mor3d, mor_r, the radius size of morphological operations on the
% binarized free-space/point-evidence (see details in include/matlab/init_freespace_clustering.m)
% pAssign, pMerge, eps_region, eps_boundary, mergin_boundary,
% mergin_region: parameters for the room segmentation (see details in include/matlab/room_segmentation.m)
% mw_min_length, alpha, noiseremoval, std: parameters for the
% preprocessiong of the input point cloud (see details in
% include/matlab/point_noise_removal.m)
scale = 1.0;
options = struct('leafsize', [13, 13, 13], 'mor3d', 0, 'mor_r', 1, 'pAssign', .9, 'pMerge', .4, 'eps_region', 10, 'eps_boundary', 10, 'mergin_boundary', 100, 'mergin_region', 100, 'mw_min_length', 800, 'alpha', 4, 'noiseremoval', 1, 'std', 1);

% Room reconstruction parameters
% See details in include/matlab/room_reconstruction_wall.m
wall_param = struct('thresh', 0, 'margin_size', 130, 'margin_sigma', 3, 'min_path_length', 0, 'penalty', 7, 'freespace_weight', 1.5, 'core_freespace_margin', 5); % margin_sigma = 3, penalty = 7 freespace = 1.5

% Wall connection analys paramters (i.e., Door addition and Room merging rules)
% See details in include/matlab/room_connection_analysis.m
door_param = struct('thresh', 600, 'wpmargin', 15, 'min_length_connectivity', 300, 'min_walldist', 1000, 'addmargin', 5);