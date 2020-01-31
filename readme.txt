Structured Indoor Modeling
Satoshi Ikehata, Hang Yan, Yasutaka Furukawa
Puslished in IEEE International Conference on Computer Vision (ICCV2015)

This package contains the matlab implementation of the partial framework in the paper above. 

Useage:
1. Compile the mex files if necessary (Manually complie each mex file or simply run include/mex/source/script_build_mex.m)
You need to install Eigen library ("http://eigen.tuxfamily.org/index.php?title=Main_Page") and put the installation path
2. Run run_apartment1.m or run_aparment2.m (you can use other dataset by changing these script and script/parameter_setup_floored.m)

Output:
- floorplan.ply (the recovered structured model in ply format); please see include/matlab/convert_floorplan_to_ply.m for more details
- floorplan.txt (the intermediate file that store the structured graph information); please see include/matlab/write_sructure_graph_basic.m for more details.
This file will be used for the future update (+detail reconstruction)

Note:
- Currently, the soure codes only recover the structured model without detail and object nodes (only recover the wall, floor, ceiling).
Other component will be added to the package later. 
- Our algorithm has a randomized process in the room segmentation step, therefore the result is not guaranteed to be same with the one in the paper.
- My code is only validated on Windows8 +　Matlab2014b

If you have questions, please contact me
Satoshi Ikeahta (sikehata@seas.wustl.edu)