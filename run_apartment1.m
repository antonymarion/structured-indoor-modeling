%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paper: Structured Indoor Modeling
% published by: IEEE International Conference on Computer Vision (ICCV) 2015
% authors: Satoshi Ikehata, Hang Yan, Yasutaka Furukawa
% contact: sikehata@wustl.edu, yanhang@seas.wustl.edu,
% furukawa@wustl.edu
% Copyright Satoshi Ikehata. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

addpath 'script';
addpath 'include/matlab';
addpath 'include/mex';

objname = 'apartment1';
datapath = 'D:/structured_modeling/structured_modeling_0.9.1/apartment1/transformed_all'; % specify the directory where all ply files are placed

print(datapath)

workspace = './workspace';
if ~exist(workspace, 'dir')
    mkdir(workspace);
end

intermediate_images = './data/apartment1/intermediate_images';
if ~exist(intermediate_images, 'dir')
    mkdir(intermediate_images);
end

dd = datestr(now,'yyyymmdd');
geometry_dir = sprintf('./data/apartment1/geometry_files/%s', dd);
if ~exist(geometry_dir, 'dir')
    mkdir(geometry_dir);
end

if ~exist(datapath, 'dir')
    error('Can not found the data directory!');
end

% Load PLY data
if ~exist(sprintf('%s/%s_data.mat', workspace, objname), 'file')
[POINT, CAM, CAMLIST] = load_ply(datapath);
save(sprintf('%s/%s_data.mat', workspace, objname), '-v7.3');
else
load(sprintf('%s/%s_data.mat', workspace, objname));
end

% Set up parameters
parameter_setup_floored
 

rng('default')

% Start structured modeling
structured_modeling




