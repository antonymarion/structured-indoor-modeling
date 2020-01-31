%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paper: Structured Indoor Modeling
% published by: IEEE International Conference on Computer Vision (ICCV) 2015
% authors: Satoshi Ikehata, Hang Yan, Yasutaka Furukawa
% contact: sikehata@wustl.edu, yanhang@seas.wustl.edu,
% furukawa@wustl.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function newevidence = normalize_evidence(evidence, alpha)
% newevidence = NORMALIZE_EVIDENCE(evidence, alpha)
% Scale the evidence from 0 to 1
% By normalizing the evidence, the evidence becomes independet of the voxel grid size

    if alpha > 0
    valid_indeces = find(evidence>0);
    values = evidence(valid_indeces);    
    mval = mean(values(:));
    stdval = std(values(:));    
    new_values = zeros(length(values),1);
    new_values(find(values>mval+alpha*stdval)) = 1;
    new_values(find(values<mval-alpha*stdval)) = 0;
    new_indeces = intersect(find(values>=mval-alpha*stdval),find(values<=mval+alpha*stdval));
    new_values(new_indeces) = 1/(2*alpha*stdval)*(values(new_indeces)-(mval-alpha*stdval));
    newevidence = zeros(size(evidence));
    newevidence(valid_indeces) = new_values;     
    else
    
    values = evidence;    
    vmax = max(values(:));
    vmin = min(values(:));
    new_values = (values-vmin)/(vmax-vmin);
    newevidence = new_values;         
    end
end

