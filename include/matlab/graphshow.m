
function graphshow(voxel_param,  WALLS, filename, doorList)
% GRAPHSHOW(WALLS, voxel_param)
% Show summarized structure graph

clustermask = zeros(voxel_param.gNum(2), voxel_param.gNum(1));
center = zeros(length(WALLS), 2);
wallcell = cell(1,1);

numRooms = length(WALLS);
for i = 1 : numRooms
    walls = WALLS{i};
    walls(end+1,:) = walls(1,:);
    x = floor(voxel_param.gNum(1)*(walls(:,1)-voxel_param.b0(1))/voxel_param.bSize(1))+1;
    y = floor(voxel_param.gNum(2)*(walls(:,2)-voxel_param.b0(2))/voxel_param.bSize(2))+1;    
    wallmask = zeros(voxel_param.gNum(2), voxel_param.gNum(1));
    for j = 1 : length(x) - 1
        x0 = min(x(j), x(j+1));
        y0 = min(y(j), y(j+1));
        x1 = max(x(j), x(j+1));
        y1 = max(y(j), y(j+1));
        xx = x0:x1;
        yy = y0:y1;
        idx = (xx-1)*voxel_param.gNum(2) + yy;
        wallmask(idx) = 1;
    end
    wallmask = imfill(wallmask,'holes');
    clustermask(find(wallmask==1)) = i;
    [a, b] = ind2sub([voxel_param.gNum(2), voxel_param.gNum(1)], find(wallmask==1));
    center(i,:) = [mean(b), mean(a)];
    
    wallcell{i} = [x, y];
%     plot(x, y, '-.w', 'LineWidth', 2)       
end

if isempty(doorList)
figure, imagesc(clustermask), axis equal, hold on

for i = 1 : numRooms
    walls = wallcell{i};
    plot(walls(:,1), walls(:,2), '-.w', 'LineWidth', 1);  
    for j = 1 : size(walls,1)-1
        wc = 0.5*(walls(j,:)+walls(j+1,:));
%         plot([center(i,1);wc(1)], [center(i,2);wc(2)], '-w')
        plot(wc(1), wc(2), 'ow', 'MarkerSize', 4, 'MarkerFaceColor', 'g');
    end    
end
plot(center(:,1), center(:,2), 'ow', 'MarkerFaceColor', 'k', 'LineWidth', 1, 'MarkerEdgeColor', 'w', 'MarkerSize', 4);
% print(gcf,'-dtiff','-r300',filename)
print('-r300','-dpng', filename);
filename2 = sprintf('%s.fig', filename);
savefig(filename2);
% saveas(gcf, filename, 'png');
else
    
    figure, imagesc(clustermask), axis equal, hold on
    WC = cell(numRooms,1);
    for i = 1 : numRooms
    walls = wallcell{i};
    plot(walls(:,1), walls(:,2), '-.w', 'LineWidth', 1);  
    for j = 1 : size(walls,1)-1
        wc(j,:) = 0.5*(walls(j,:)+walls(j+1,:));
%         plot([center(i,1);wc(j,1)], [center(i,2);wc(j,2)], '-w')
%         plot(wc(1), wc(2), 'ow', 'MarkerSize', 4, 'MarkerFaceColor', 'g');
    end    
    WC{i} = wc;
    end
    

        
    for i = 1 : numRooms
        wc = WC{i};
    for j = 1 : size(wc,1)
        plot(wc(j,1), wc(j,2), 'ow', 'MarkerSize', 4, 'MarkerFaceColor', 'g');
    end    
    end     
    plot(center(:,1), center(:,2), 'ow', 'MarkerFaceColor', 'k', 'LineWidth', 1, 'MarkerEdgeColor', 'w', 'MarkerSize', 4);
    
        doorcenter = zeros(length(doorList),2);
    for k = 1 : length(doorList)
        if doorList{k}.roomId1 == 0 || doorList{k}.roomId2 == 0
            return;
        end
        if doorList{k}.aisleFlag == 0

            x0 = WC{doorList{k}.roomId1}(doorList{k}.wallId1,:);
            x1 = WC{doorList{k}.roomId2}(doorList{k}.wallId2,:);
            doorcenter(k,:) = 0.5*(x0+x1);
            plot([x0(1);doorcenter(k,1)], [x0(2);doorcenter(k,2)], '-r');
            plot([x1(1);doorcenter(k,1)], [x1(2);doorcenter(k,2)], '-r');            
        else
            x0 = WC{doorList{k}.roomId1}(doorList{k}.wallId1,:);
            x1 = WC{doorList{k}.roomId2}(doorList{k}.wallId2,:);
            plot([x0(1);x1(1)], [x0(2);x1(2)], '-k');
        end        
    end
    plot(doorcenter(:,1), doorcenter(:,2), 'sw', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
   
    print('-r300','-dpng', filename);
    filename2 = sprintf('%s.fig', filename);
    savefig(filename2);
    
end

end

