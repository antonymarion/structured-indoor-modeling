function [labels connectivity] = find_connectivity(cluster_fs)
[h, w] = size(cluster_fs);
labels = unique(cluster_fs);
labels = labels(labels~=0);
connectivity = cell(length(labels),1);
for i = 1 : length(labels)
    B = cell2mat(bwboundaries(cluster_fs==labels(i)));
    for j = 1 : size(B,1)
        y = B(j,1);
        x = B(j,2);
        if(x+1 <= w)
        if cluster_fs(y,x+1) ~= labels(i) & cluster_fs(y,x+1) ~= 0
            connectivity{i,1}(end+1) =  cluster_fs(y,x+1);      
        end
        end
        if(x-1 > 0)
        if cluster_fs(y,x-1) ~= labels(i) & cluster_fs(y,x-1) ~= 0
            connectivity{i,1}(end+1) =  cluster_fs(y,x-1);      
        end
        end
        if(y+1 <= h)
        if cluster_fs(y+1,x) ~= labels(i) & cluster_fs(y+1,x) ~= 0
            connectivity{i,1}(end+1) =  cluster_fs(y+1,x);      
        end
        end
        
        if(y-1>0)
        if cluster_fs(y-1,x) ~= labels(i) & cluster_fs(y-1,x) ~= 0
            connectivity{i,1}(end+1) =  cluster_fs(y-1,x);            
        end
        end
    end
    connectivity{i} = unique(connectivity{i});
end