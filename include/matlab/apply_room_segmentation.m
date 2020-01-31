%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paper: Structured Indoor Modeling
% published by: IEEE International Conference on Computer Vision (ICCV) 2015
% authors: Satoshi Ikehata, Hang Yan, Yasutaka Furukawa
% contact: sikehata@wustl.edu, yanhang@seas.wustl.edu,
% furukawa@wustl.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [minIndex, minDistance, ambiguity, clusterCenters, indexOnBoundary] = apply_room_segmentation(feature, BW, init_cluster_indeces, indexRegion, maxIter, pAssign, pMerge, xBoundary, yBoundary)
% [minIndex, minDistance, ambiguity, clusterCenters, indexOnBoundary] = APPLY_ROOM_SEGMENTATION(feature, BW, init_cluster_indeces, indexRegion, maxIter, pAssign, pMerge, xBoundary, yBoundary)
% Room segmentation by K-medoids algorithm

h = figure;
indexOnBoundary = zeros(length(xBoundary),1);
[numRegion, numBoundary] = size(feature);

indexCentersOnRegion = init_cluster_indeces;

% initialize cluster centers
centers = feature(init_cluster_indeces,:);
indexCenters = indexRegion(init_cluster_indeces);

% for stopping criteria
ListNumClusters = zeros(3,1);

% main step of k-medoids clustering
numClusters = size(centers,1);
ListNumClusters(1) = numClusters;
ListNumClustersIndex = 1;
for iter = 1 : maxIter
    
display(sprintf('Iter %d:  K-medoids clustering (#cluster is %d)',iter,numClusters))
    
%% assign new cluster

minDistance = 1.0e12*ones(numRegion,1);
Distance2nd = 1.0e12*ones(numRegion,1);
minIndex = ones(numRegion,1);
sum_target = sum(feature');
sum_target = sum_target(:);

for k = 1 : numClusters   
    sum_seed = sum(centers(k,:)');
    sum_seed = sum_seed(:);      
    hammingDistance = mex_dist_weighted_double((feature'),(centers(k,:)'),sum_target, sum_seed)/2;
    % if all elements in a featue are zero, the maximum hamming
    hammingDistance(sum_seed==0) = 1;
    hammingDistance(sum_target==0) = 1;   
    id = find(hammingDistance < minDistance);
    id2 = find(hammingDistance >= minDistance & hammingDistance < Distance2nd);
    Distance2nd(id2) = hammingDistance(id2);
    Distance2nd(id) = minDistance(id);
    minDistance(id) = hammingDistance(id);   
    minIndex(id) = k;   
end
ambiguity =  minDistance./Distance2nd;
ambiguity = ambiguity/max(ambiguity(:));        

[h, w] = size(BW);  

% delete empty clusters
isValid = logical(ones(numClusters,1));
cnt = 0;
for k = 1 : numClusters
    indexKthCluster = find(minIndex == k);
    if isempty(indexKthCluster)
       isValid(k) = 0;
       cnt = cnt + 1;
    else
       minIndex(find(minIndex == k)) = k - cnt;            
    end
end
centers = centers(isValid,:);
indexCenters = indexCenters(isValid);
indexCentersOnRegion =  indexCentersOnRegion(isValid);
numClusters = size(centers,1);
isValid = logical(ones(numClusters,1));

%% update cluster center    
%     % update cluster center (exact)
%     display('...Updating centers')
%     new_centers = zeros(size(centers));
%     for k = 1 : numClusters
%         indexKthCluster = find(minIndex == k);
%         cluster_features = feature(indexKthCluster,:);
%         sum_target = sum(logical(cluster_features'));
%         sum_target = sum_target(:);
%         sum_seed = sum(logical(cluster_features'));
%         sum_seed = sum_seed(:);
%         hammingDistance = mex_dist_weighted_hamming(logical(cluster_features'),logical(cluster_features'),sum_target, sum_seed)/2;
%         hammingDistance(sum_seed==0) = 1;
%         hammingDistance(sum_target==0) = 1;         
%         hammingDistance = reshape(hammingDistance, size(cluster_features,1), size(cluster_features,1));
%         [~, newCenterIndex] = min(sum(hammingDistance));
%         new_centers(k,:) = cluster_features(newCenterIndex,:);
%         indexCenters(k) = indexRegion(indexKthCluster(newCenterIndex));
%     end

% update cluster center (approximate) fast, but it is not guaranteed that the potential is
%     converged into the global minima
display('...Updating centers')
new_centers = zeros(size(centers));
for k = 1 : numClusters
    indexKthCluster = find(minIndex == k);
    cluster_features = feature(indexKthCluster,:);

    if size(cluster_features,1) > 40 + numClusters            
        currentCost = 1.0e12;
        currentCenterIndex = -1;
        sum_target = sum((cluster_features'));
        sum_target = sum_target(:);
        for i = 1 : 10
            % only extract (40+k) samples
            data_num=[1:size(cluster_features,1)];
            ind = randperm(length(data_num));
            subIndex = ind(1:40+numClusters);
            subClusterFeatures = cluster_features(subIndex,:);
            sum_target = sum((subClusterFeatures'));
            sum_target = sum_target(:);
            sum_seed = sum((subClusterFeatures'));
            sum_seed = sum_seed(:);
            hammingDistance = mex_dist_weighted_double((subClusterFeatures'),(subClusterFeatures'),sum_target, sum_seed)/2;
            hammingDistance(sum_seed==0) = 1;
            hammingDistance(sum_target==0) = 1;         
            hammingDistance = reshape(hammingDistance, size(subClusterFeatures,1), size(subClusterFeatures,1));
            [~, subCenterIndex] = min(sum(hammingDistance));
            newCenterIndex = subIndex(subCenterIndex);               

            sum_seed = sum((cluster_features(newCenterIndex,:)'));
            sum_seed = sum_seed(:);
            hammingDistance = mex_dist_weighted_double((subClusterFeatures'),(cluster_features(newCenterIndex,:)'),sum_target, sum_seed)/2;           
            hammingDistance(sum_seed==0) = 1;
            hammingDistance(sum_target==0) = 1;         
            newCost = sum(hammingDistance);         
            if newCost < currentCost
                currentCenterIndex = newCenterIndex;
                currentCost = newCost;
            end                
        end

        newCenterIndex = currentCenterIndex;
        new_centers(k,:) = cluster_features(newCenterIndex,:);
        indexCenters(k) = indexRegion(indexKthCluster(newCenterIndex)); 
        indexCentersOnRegion(k) = indexKthCluster(newCenterIndex);

    else         
    sum_target = sum((cluster_features'));
    sum_target = sum_target(:);
    sum_seed = sum((cluster_features'));
    sum_seed = sum_seed(:);
    hammingDistance = mex_dist_weighted_double((cluster_features'),(cluster_features'),sum_target, sum_seed)/2;
    hammingDistance(sum_seed==0) = 1;
    hammingDistance(sum_target==0) = 1;         
    hammingDistance = reshape(hammingDistance, size(cluster_features,1), size(cluster_features,1));
    [~, newCenterIndex] = min(sum(hammingDistance));
    new_centers(k,:) = cluster_features(newCenterIndex,:);
    indexCenters(k) = indexRegion(indexKthCluster(newCenterIndex));
    indexCentersOnRegion(k) = indexKthCluster(newCenterIndex);
    end
end

% visualization
color = jet(2*numClusters);
data_num=[1:2*numClusters];
ind = randperm(length(data_num));
color = color(ind,:);
color = color(1:numClusters,:);

[h, w] = size(BW);

figure(h)
tempResult = zeros(h,w);
tempResult(BW>0) = 1;

imagesc(tempResult), axis equal, hold on
uniqueindex = unique(minIndex);
for tt = 1 : length(uniqueindex)
    [yy,xx] = ind2sub([h, w], indexRegion(find(minIndex==uniqueindex(tt))));    
    plot(xx, yy, '.', 'MarkerEdgeColor', color(uniqueindex(tt),:), 'MarkerSize', 5);
end

[y, x] = ind2sub([h,w],indexCenters);
scatter(x, y, 48, color(1:numClusters,:),'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w');    



% stopping criteria
ListNumClusters(mod(ListNumClustersIndex,3)+1) = numClusters;
ListNumClustersIndex = ListNumClustersIndex + 1;

if ListNumClusters(1) ==  ListNumClusters(2) & ListNumClusters(1) ==  ListNumClusters(3) & iter > 5
   break;        
end

if numClusters == 1
   break; 
end
centers = new_centers;      

%% merge clusters & upsate centers again
% only merge nearest two clusters
display('...Merging clusters');    
count = 1;
while (1)     
sum_target = sum((centers)');
sum_target = sum_target(:);
sum_seed = sum((centers)');
sum_seed = sum_seed(:);
hammingDistance = mex_dist_weighted_double((centers)',(centers)',sum_target, sum_seed)/2;
hammingDistance(sum_seed==0) = 1;
hammingDistance(sum_target==0) = 1;
hammingDistance = reshape(hammingDistance, numClusters, numClusters);

[x y] = meshgrid(1:numClusters, 1:numClusters);
isNonDiagonal = find(x ~= y);
[minDist, minID] = min(hammingDistance(isNonDiagonal));
mergeIndexPair = isNonDiagonal(minID);  
if minDist > pMerge
    break;
end

if isempty(mergeIndexPair)
    break;
end
[mergeIndex1 mergeIndex2] = ind2sub([numClusters numClusters], mergeIndexPair);
if mergeIndex1 > mergeIndex2
    temp = mergeIndex1;
    mergeIndex1 = mergeIndex2;
    mergeIndex2 = temp;
end  

% center updates 
isValid(mergeIndex2) = 0;

indexPairClusters = [find(minIndex==mergeIndex1);find(minIndex == mergeIndex2)];         
cluster_features = feature(indexPairClusters,:);


%     % exact
%     sum_target = sum(logical(cluster_features'));
%     sum_target = sum_target(:);
%     sum_seed = sum(logical(cluster_features'));
%     sum_seed = sum_seed(:);
%     hammingDistance = mex_dist_weighted_hamming(logical(cluster_features'),logical(cluster_features'), sum_target, sum_seed)/2;
%     hammingDistance(sum_seed==0) = 1;
%     hammingDistance(sum_target==0) = 1;
%     hammingDistance = reshape(hammingDistance, size(cluster_features,1), size(cluster_features,1));
%     [minDistance, newCenterIndex] = min(sum(hammingDistance));
%     centers(mergeIndex1,:) = cluster_features(newCenterIndex,:);
%     indexCenters(mergeIndex1) = indexRegion(indexPairClusters(newCenterIndex));  
%     centers = centers(isValid,:);    
%     indexCenters = indexCenters(isValid);
%     numClusters = size(centers,1);    

% approximate (randomized sarch, motivated from CLARANS)
current_center = centers(mergeIndex1,:);
current_id = -1;
sum_target = sum((cluster_features'));
sum_target = sum_target(:);
sum_seed = sum((current_center'));
sum_seed = sum_seed(:);
hammingDistance = mex_dist_weighted_double((cluster_features'),(current_center'), sum_target, sum_seed)/2;
hammingDistance(sum_seed==0) = 1;
hammingDistance(sum_target==0) = 1;

data_num=[1:size(cluster_features,1)];
maxNeighbor = floor((size(cluster_features,1)-1)/80);

current_cost = size(cluster_features,1);
ind = randperm(length(data_num));    
for i = 1 : maxNeighbor
    new_center = cluster_features(ind(i),:);
    sum_target = sum((cluster_features'));
    sum_target = sum_target(:);
    sum_seed = sum((new_center'));
    sum_seed = sum_seed(:);
    hammingDistance = mex_dist_weighted_double((cluster_features'),(new_center'), sum_target, sum_seed)/2;
    hammingDistance(sum_seed==0) = 1;
    hammingDistance(sum_target==0) = 1;
    new_cost = sum(hammingDistance(:));  
    if new_cost < current_cost
        current_center = new_center;
        current_cost = new_cost;
        current_id = ind(i);
    end
end

centers(mergeIndex1,:) = current_center;
if current_id ~= -1
indexCenters(mergeIndex1) = indexRegion(indexPairClusters(current_id)); 
indexCentersOnRegion(mergeIndex1) = indexPairClusters(current_id);
end
centers = centers(isValid,:);    
indexCenters = indexCenters(isValid);
indexCentersOnRegion =  indexCentersOnRegion(isValid);
numClusters = size(centers,1);     

%   
minIndex(find(minIndex == mergeIndex2)) = mergeIndex1;
minIndex(find(minIndex > mergeIndex2)) = minIndex(find(minIndex > mergeIndex2)) - 1;
isValid = logical(ones(numClusters,1));       
count = count + 1;
end
[h, w] = size(BW);
[y, x] = ind2sub([h,w],indexCenters);
clusterCenters = [x,y];

% visualization
color = jet(numClusters);
data_num=[1:numClusters];
ind = randperm(length(data_num));
color = color(ind,:);

[h, w] = size(BW);
figure(h)
tempResult = zeros(h*w,3);
tempResult(:,1) = BW(:);
tempResult(:,2) = BW(:);
tempResult(:,3) = BW(:);
tempResult(indexRegion,1) = color(minIndex,1);
tempResult(indexRegion,2) = color(minIndex,2);
tempResult(indexRegion,3) = color(minIndex,3);
imagesc(reshape(tempResult,h, w, 3))
axis equal
hold on
[y, x] = ind2sub([h,w],indexCenters);
scatter(x, y, 64, color,'MarkerFaceColor', 'k');    

pause(0.5)  
display(sprintf('K-medios distance = %.2f', sum(hammingDistance(:))));
  
end
