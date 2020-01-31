function newmap = deleteIsolatedCluster(map)

[h w] = size(map);
newmap = zeros(h, w); 
clusterIndex = unique(map);
count = 1;
for i = 1 : length(clusterIndex)
if(clusterIndex(i) ~= 0)
indexKthCluster = find(map==clusterIndex(i));
binaryCluster = zeros(h, w);
binaryCluster(indexKthCluster) = 1;
CC = bwconncomp(binaryCluster,4);
numPixels = cellfun(@numel,CC.PixelIdxList);
[biggest,idx] = max(numPixels);
newmap(CC.PixelIdxList{idx}) = count;
count = count + 1;
end
end

end