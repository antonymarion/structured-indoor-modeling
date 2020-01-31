%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paper: Structured Indoor Modeling
% published by: IEEE International Conference on Computer Vision (ICCV) 2015
% authors: Satoshi Ikehata, Hang Yan, Yasutaka Furukawa
% contact: sikehata@wustl.edu, yanhang@seas.wustl.edu,
% furukawa@wustl.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [doorList] = check_wall_connection_type(connectivity, freespace_evidence3d, point_evidence3d, voxel_param, options, door_param, WALLS, LABELS)
% [doorList] = CHECK_WALL_CONNECTION_TYPE(connectivity, freespace_evidence3d, point_evidence3d, voxel_param, options, door_param, WALLS, LABELS);
% Pick a pair of walls from different room nodes and then check if any
% connections
% For efficiencty and robustness, we implemented to check a pair of walls that satisfy follows
% (i) parent room nodes share the boundary in the segmennted 2-D free-space evidence
% (ii) walls are facing with each other
% (iii) the spatial distance of two walls are less than thre shold

close all

doorList = cell(1,1);
numRooms = size(WALLS,1);

freespace_evidence = reshape(freespace_evidence3d, voxel_param.gNum(2)*voxel_param.gNum(1), voxel_param.gNum(3));
f3d = freespace_evidence;
point_evidence = reshape(point_evidence3d, voxel_param.gNum(1)*voxel_param.gNum(2), voxel_param.gNum(3));
p3d = point_evidence;
point_evidence = sum((point_evidence>0)');
point_evidence = reshape(point_evidence, voxel_param.gNum(2), voxel_param.gNum(1));
point_evidence = normalize_evidence(point_evidence, 1);

%% room pair candidate

doorCount = 1;
connectivityMap = zeros(numRooms, numRooms);
for i = 1 : numRooms
    numConnectivity = length(connectivity{i});
    for k = 1 : numConnectivity
        j = connectivity{i}(k);
        if connectivityMap(i,j) == 0
       
        connectivityMap(i,j) = 1;
        connectivityMap(j,i) = 1;
        
        wall0 = WALLS{i};
        wall0(end+1,:) = wall0(1,:);
        wall1 = WALLS{j};
        wall1(end+1,:) = wall1(1,:);
        label0 = LABELS{i};
        label1 = LABELS{j};        
        % for visualization
        x0 = floor(voxel_param.gNum(1)*(wall0(:,1)-voxel_param.b0(1))/voxel_param.bSize(1))+1;
        y0 = floor(voxel_param.gNum(2)*(wall0(:,2)-voxel_param.b0(2))/voxel_param.bSize(2))+1;
        x0(end+1) = x0(1);
        y0(end+1) = y0(1); 
        x1 = floor(voxel_param.gNum(1)*(wall1(:,1)-voxel_param.b0(1))/voxel_param.bSize(1))+1;
        y1 = floor(voxel_param.gNum(2)*(wall1(:,2)-voxel_param.b0(2))/voxel_param.bSize(2))+1;
        x1(end+1) = x1(1);
        y1(end+1) = y1(1); 
        

        
        for p = 1 : size(wall0,1)-1
            for q = 1 : size(wall1,1)-1

              % x direction
              % first top second bottom
                if label0(p) == 1 & label1(q) == -1 & max(wall0(p,2), wall0(p+1,2)) < min(wall1(q,2), wall1(q+1,2))
                    pw1 = floor(mean([y0(p) y0(p+1)]));
                    pw2 = floor(mean([y1(q) y1(q+1)]));
                    w = min(door_param.wpmargin, abs(floor(0.5*(pw1-pw2))));   

                    % first
                    wallProfile0 = zeros(voxel_param.gNum(3), voxel_param.gNum(1));
                    freeSpaceProfile0 = ones(voxel_param.gNum(3), voxel_param.gNum(1));
                    xmin = floor(min(x0(p), x0(p+1)))-door_param.addmargin;
                    xmax = floor(max(x0(p), x0(p+1)))+door_param.addmargin;
                    xmin = max(xmin,1);
                    xmax = min(xmax,voxel_param.gNum(1));                    
                    yy = floor(mean([y0(p) y0(p+1)]));

                    for tt = -door_param.wpmargin : w                                                
                          t = [xmin:xmax];
                          [Z, X] = meshgrid([1:voxel_param.gNum(3)], t);                 
                          id = (Z-1).*voxel_param.gNum(2)*voxel_param.gNum(1) + (X-1)*voxel_param.gNum(2) + yy + tt;
                          f2d = zeros(size(freeSpaceProfile0));
                          f2d((X-1).*voxel_param.gNum(3)+Z) = f3d(id);
                          f2d = imclose(f2d, strel('disk', 5)); 

                          p2d = zeros(size(wallProfile0));
                          p2d((X-1).*voxel_param.gNum(3)+Z) = p3d(id);
                          f2d = imclose(f2d, strel('disk', 5));  

                          freeSpaceProfile0 = min(freeSpaceProfile0, f2d);
                          wallProfile0 =  wallProfile0 + p2d;
                    end 
                    wallProfile0 = normalize_evidence(wallProfile0, 0);
                    freeSpaceProfile0 = normalize_evidence(freeSpaceProfile0, 0);


                    % second
                    wallProfile1 = zeros(voxel_param.gNum(3), voxel_param.gNum(1));
                    freeSpaceProfile1 = ones(voxel_param.gNum(3), voxel_param.gNum(1));
                    xmin = floor(min(x1(q), x1(q+1)))-door_param.addmargin;
                    xmax = floor(max(x1(q), x1(q+1)))+door_param.addmargin;
                    xmin = max(xmin,1);
                    xmax = min(xmax,voxel_param.gNum(1));                    
                    yy = floor(mean([y1(q) y1(q+1)]));

                    for tt = -w : door_param.wpmargin                                                
                          t = [xmin:xmax];
                          [Z, X] = meshgrid([1:voxel_param.gNum(3)], t);                 
                          id = (Z-1).*voxel_param.gNum(2)*voxel_param.gNum(1) + (X-1)*voxel_param.gNum(2) + yy + tt;
                          f2d = zeros(size(freeSpaceProfile1));
                          f2d((X-1).*voxel_param.gNum(3)+Z) = f3d(id);
                          f2d = imclose(f2d, strel('disk', 5)); 

                          p2d = zeros(size(wallProfile1));
                          p2d((X-1).*voxel_param.gNum(3)+Z) = p3d(id);
                          f2d = imclose(f2d, strel('disk', 5));  

                          freeSpaceProfile1 = min(freeSpaceProfile1, f2d);
                          wallProfile1 =  wallProfile1 + p2d;
                    end 
                    wallProfile1 = normalize_evidence(wallProfile1, 0);
                    freeSpaceProfile1 = normalize_evidence(freeSpaceProfile1, 0);
                    
                    if abs(pw1-pw2) < door_param.min_walldist/options.leafsize(1)
                    min_x = max(min(x0(p), x0(p+1)),  min(x1(q), x1(q+1)));
                    max_x = min(max(x0(p), x0(p+1)),  max(x1(q), x1(q+1)));
  
                    [doorX0, doorX1, minCost, flag] = wall_profile_analysis(freeSpaceProfile0, freeSpaceProfile1, wallProfile0, wallProfile1, min_x, max_x, door_param.thresh/options.leafsize(1), i, j, p, q); % doorX, doorY: c++coordinate
                    if flag ~= 0
                        mind = max(abs(doorX0(1) - min_x), abs(doorX1(1) - max_x));
                        figure('WindowStyle', 'docked', 'Name', sprintf('Door is between [Room %d Wall %d] [Room %d Wall %d]', i,p,j,q)),
                        imagesc(point_evidence), axis equal
                        hold on
                        plot([x0(p);x0(p+1)],[y0(p);y0(p+1)],'o-r', 'LineWidth',3)     
                        plot([x1(q);x1(q+1)],[y1(q);y1(q+1)],'o-g', 'LineWidth',3)    
                        title(sprintf('%f %d %d %f %f %f %f', mind, label0(p), label1(q), wall0(p,1), wall0(p,2), wall1(q,1), wall1(q,2)))
                        hold off
                        axis equal

                        xDoor0 = voxel_param.bSize(1)*doorX0(1)/voxel_param.gNum(1)+voxel_param.b0(1);
                        zDoor0 = voxel_param.bSize(3)*doorX0(2)/voxel_param.gNum(3)+voxel_param.b0(3); 
                        xDoor1 = voxel_param.bSize(1)*doorX1(1)/voxel_param.gNum(1)+voxel_param.b0(1);
                        zDoor1 = voxel_param.bSize(3)*doorX1(2)/voxel_param.gNum(3)+voxel_param.b0(3); 
                        yDoor0 = mean([wall0(p,2), wall0(p+1,2)]); 
                        yDoor1 = mean([wall1(q,2), wall1(q+1,2)]);                                                                 
                        
                        aisleFlag = 0;
                        if flag == 2
                           aisleFlag = 1; 
                        end
                        if flag > 0
                        doorINFO = struct('aisleFlag', aisleFlag, 'x0', [xDoor0, -1, zDoor0], 'x1', [xDoor1, -1, zDoor1], 'roomId1', i, 'roomId2', j, 'wallId1', p, 'wallId2', q);
                        doorList{doorCount} = doorINFO;
                        doorCount = doorCount + 1;     
                        end
                    end 
                    end
                end
              
                % x dir
                % first bottom second top
                if label0(p) == -1 & label1(q) == 1  & min(wall0(p,2), wall0(p+1,2)) > max(wall1(q,2), wall1(q+1,2))                 
                    pw1 = floor(mean([y0(p) y0(p+1)]));
                    pw2 = floor(mean([y1(q) y1(q+1)]));
                    w = min(door_param.wpmargin, abs(floor(0.5*(pw1-pw2))));   

                    % first
                    wallProfile0 = zeros(voxel_param.gNum(3), voxel_param.gNum(1));
                    freeSpaceProfile0 = ones(voxel_param.gNum(3), voxel_param.gNum(1));
                    xmin = floor(min(x0(p), x0(p+1)))-door_param.addmargin;
                    xmax = floor(max(x0(p), x0(p+1)))+door_param.addmargin;
                    xmin = max(xmin,1);
                    xmax = min(xmax,voxel_param.gNum(1));                    
                    yy = floor(mean([y0(p) y0(p+1)]));

                    for tt = -w : door_param.wpmargin                                                
                          t = [xmin:xmax];
                          [Z, X] = meshgrid([1:voxel_param.gNum(3)], t);                 
                          id = (Z-1).*voxel_param.gNum(2)*voxel_param.gNum(1) + (X-1)*voxel_param.gNum(2) + yy + tt;
                          f2d = zeros(size(freeSpaceProfile0));
                          f2d((X-1).*voxel_param.gNum(3)+Z) = f3d(id);
                          f2d = imclose(f2d, strel('disk', 5)); 

                          p2d = zeros(size(wallProfile0));
                          p2d((X-1).*voxel_param.gNum(3)+Z) = p3d(id);
                          f2d = imclose(f2d, strel('disk', 5));  

                          freeSpaceProfile0 = min(freeSpaceProfile0, f2d);
                          wallProfile0 =  wallProfile0 + p2d;
                    end 
                    wallProfile0 = normalize_evidence(wallProfile0, 0);
                    freeSpaceProfile0 = normalize_evidence(freeSpaceProfile0, 0);


                    % second
                    wallProfile1 = zeros(voxel_param.gNum(3), voxel_param.gNum(1));
                    freeSpaceProfile1 = ones(voxel_param.gNum(3), voxel_param.gNum(1));
                    xmin = floor(min(x1(q), x1(q+1)))-door_param.addmargin;
                    xmax = floor(max(x1(q), x1(q+1)))+door_param.addmargin;
                    xmin = max(xmin,1);
                    xmax = min(xmax,voxel_param.gNum(1));                    
                    yy = floor(mean([y1(q) y1(q+1)]));

                    for tt = -door_param.wpmargin : w                                                
                          t = [xmin:xmax];
                          [Z, X] = meshgrid([1:voxel_param.gNum(3)], t);                 
                          id = (Z-1).*voxel_param.gNum(2)*voxel_param.gNum(1) + (X-1)*voxel_param.gNum(2) + yy + tt;
                          f2d = zeros(size(freeSpaceProfile1));
                          f2d((X-1).*voxel_param.gNum(3)+Z) = f3d(id);
                          f2d = imclose(f2d, strel('disk', 5)); 

                          p2d = zeros(size(wallProfile1));
                          p2d((X-1).*voxel_param.gNum(3)+Z) = p3d(id);
                          f2d = imclose(f2d, strel('disk', 5));  

                          freeSpaceProfile1 = min(freeSpaceProfile1, f2d);
                          wallProfile1 =  wallProfile1 + p2d;
                    end 
                    wallProfile1 = normalize_evidence(wallProfile1, 0);
                    freeSpaceProfile1 = normalize_evidence(freeSpaceProfile1, 0);
                    
                    
                    if abs(pw1-pw2) < door_param.min_walldist/options.leafsize(1)
                    min_x = max(min(x0(p), x0(p+1)),  min(x1(q), x1(q+1)));
                    max_x = min(max(x0(p), x0(p+1)),  max(x1(q), x1(q+1)));
                    [doorX0, doorX1, minCost, flag] = wall_profile_analysis(freeSpaceProfile0, freeSpaceProfile1, wallProfile0, wallProfile1, min_x, max_x, door_param.thresh/options.leafsize(1), i, j, p, q);
                     if flag ~= 0                  
                         mind = max(abs(doorX0(1) - min_x), abs(doorX1(1) - max_x));
                        figure('WindowStyle', 'docked', 'Name', sprintf('Door is between [Room %d Wall %d] [Room %d Wall %d]', i,p,j,q)),
                        imagesc(point_evidence), axis equal
                        hold on
                        plot([x0(p);x0(p+1)],[y0(p);y0(p+1)],'o-r', 'LineWidth',3)     
                        plot([x1(q);x1(q+1)],[y1(q);y1(q+1)],'o-g', 'LineWidth',3)     
                        title(sprintf('%f %d %d %f %f %f %f', mind, label0(p), label1(q), wall0(p,1), wall0(p,2), wall1(q,1), wall1(q,2)))
                        hold off
                        axis equal
                        
                        xDoor0 = voxel_param.bSize(1)*doorX0(1)/voxel_param.gNum(1)+voxel_param.b0(1);
                        zDoor0 = voxel_param.bSize(3)*doorX0(2)/voxel_param.gNum(3)+voxel_param.b0(3); 
                        xDoor1 = voxel_param.bSize(1)*doorX1(1)/voxel_param.gNum(1)+voxel_param.b0(1);
                        zDoor1 = voxel_param.bSize(3)*doorX1(2)/voxel_param.gNum(3)+voxel_param.b0(3); 
                        yDoor0 = mean([wall0(p,2), wall0(p+1,2)]); 
                        yDoor1 = mean([wall1(q,2), wall1(q+1,2)]);      
                    
                        aisleFlag = 0;
                        
                        if flag == 2
                           aisleFlag = 1; 
                        end
                        if flag > 0
                        doorINFO = struct('aisleFlag', aisleFlag, 'x0', [xDoor0, -1, zDoor0], 'x1', [xDoor1, -1, zDoor1], 'roomId1', i, 'roomId2', j, 'wallId1', p, 'wallId2', q);
                        doorList{doorCount} = doorINFO;
                        doorCount = doorCount + 1;  
                        end
                    end                                      
                    end
                end
             
                % ydir
                % first right second left
                if label0(p) == 2 & label1(q) == -2 & max(wall1(q,1), wall1(q+1,1)) < min(wall0(p,1),wall0(p+1,1)) 
                   
                    pw1 = floor(mean([x0(p) x0(p+1)]));
                    pw2 = floor(mean([x1(q) x1(q+1)]));
                    w = min(door_param.wpmargin, abs(floor(0.5*(pw1-pw2))));   
                    
                    % first
                    wallProfile0 = zeros(voxel_param.gNum(3), voxel_param.gNum(2));
                    freeSpaceProfile0 = ones(voxel_param.gNum(3), voxel_param.gNum(2));
                    ymin = floor(min(y0(p), y0(p+1)))-door_param.addmargin;
                    ymax = floor(max(y0(p), y0(p+1)))+door_param.addmargin;
                    ymin = max(ymin,1);
                    ymax = min(ymax,voxel_param.gNum(2));                    
                    xx = floor(mean([x0(p) x0(p+1)]));
                    
                    for tt = -w : door_param.wpmargin                                                
                          t = [ymin:ymax];
                          [Z, Y] = meshgrid([1:voxel_param.gNum(3)], t);                 
                          id = (Z-1).*voxel_param.gNum(2)*voxel_param.gNum(1) + (xx+tt-1)*voxel_param.gNum(2) + Y;
                          f2d = zeros(size(freeSpaceProfile0));
                          f2d((Y-1).*voxel_param.gNum(3)+Z) = f3d(id);
                          f2d = imclose(f2d, strel('disk', 5)); 
                          
                          p2d = zeros(size(wallProfile0));
                          p2d((Y-1).*voxel_param.gNum(3)+Z) = p3d(id);
                          
                          f2d = imclose(f2d, strel('disk', 5)); 
%                           p2d = imclose(p2d, strel('disk', 5)); 
                          
                          freeSpaceProfile0 = min(freeSpaceProfile0, f2d);
                          wallProfile0 =  wallProfile0 + p2d;
                    end 
                    wallProfile0 = normalize_evidence(wallProfile0, 0);
                    freeSpaceProfile0 = normalize_evidence(freeSpaceProfile0, 0);
                    
                    
                    % second
                    wallProfile1 = zeros(voxel_param.gNum(3), voxel_param.gNum(2));
                    freeSpaceProfile1 = ones(voxel_param.gNum(3), voxel_param.gNum(2));
                    ymin = floor(min(y1(q), y1(q+1)))-door_param.addmargin;
                    ymax = floor(max(y1(q), y1(q+1)))+door_param.addmargin;
                    ymin = max(ymin,1);
                    ymax = min(ymax,voxel_param.gNum(2));                    
                    xx = floor(mean([x1(q) x1(q+1)]));
                    
                    for tt = -door_param.wpmargin : w                                                
                          t = [ymin:ymax];
                          [Z, Y] = meshgrid([1:voxel_param.gNum(3)], t);                 
                          id = (Z-1).*voxel_param.gNum(2)*voxel_param.gNum(1) + (xx+tt-1)*voxel_param.gNum(2) + Y;
                          f2d = zeros(size(freeSpaceProfile1));
                          f2d((Y-1).*voxel_param.gNum(3)+Z) = f3d(id);

                          p2d = zeros(size(wallProfile1));
                          p2d((Y-1).*voxel_param.gNum(3)+Z) = p3d(id);                          
                          f2d = imclose(f2d, strel('disk', 5));                           
                          freeSpaceProfile1 = min(freeSpaceProfile1, f2d);
                          wallProfile1 =  wallProfile1 + p2d;
                    end 
                    wallProfile1 = normalize_evidence(wallProfile1, 0);
                    freeSpaceProfile1 = normalize_evidence(freeSpaceProfile1, 0);  
                    
                    if abs(pw1-pw2) < door_param.min_walldist/options.leafsize(1)
                    min_x = max(min(y0(p), y0(p+1)),  min(y1(q), y1(q+1)));
                    max_x = min(max(y0(p), y0(p+1)),  max(y1(q), y1(q+1)));
                    [doorY0, doorY1, minCost, flag] = wall_profile_analysis(freeSpaceProfile0, freeSpaceProfile1, wallProfile0, wallProfile1, min_x, max_x, door_param.thresh/options.leafsize(1), i, j, p, q);

                     if flag ~= 0
                        mind = max(abs(doorY0(1) - min_x), abs(doorY1(1) - max_x));
                        figure('WindowStyle', 'docked', 'Name', sprintf('Door is between [Room %d Wall %d] [Room %d Wall %d]', i,p,j,q)),
                        imagesc(point_evidence), axis equal
                        hold on
                        plot([x0(p);x0(p+1)],[y0(p);y0(p+1)],'o-r', 'LineWidth',3)     
                        plot([x1(q);x1(q+1)],[y1(q);y1(q+1)],'o-g', 'LineWidth',3)       
                        title(sprintf('%f %d %d %f %f %f %f', mind, label0(p), label1(q), wall0(p,1), wall0(p,2), wall1(q,1), wall1(q,2)))
                        hold off
                        axis equal
                        
                        yDoor0 = voxel_param.bSize(2)*doorY0(1)/voxel_param.gNum(2)+voxel_param.b0(2);
                        zDoor0 = voxel_param.bSize(3)*doorY0(2)/voxel_param.gNum(3)+voxel_param.b0(3); 
                        yDoor1 = voxel_param.bSize(2)*doorY1(1)/voxel_param.gNum(2)+voxel_param.b0(2);
                        zDoor1 = voxel_param.bSize(3)*doorY1(2)/voxel_param.gNum(3)+voxel_param.b0(3); 
                        xDoor0 = mean([wall0(p,1), wall0(p+1,1)]); 
                        xDoor1 = mean([wall1(q,1), wall1(q+1,1)]);                                   
                   
                        
                        aisleFlag = 0;
                        if flag == 2
                           aisleFlag = 1; 
                        end
                        if flag > 0
                        doorINFO = struct('aisleFlag', aisleFlag, 'x0', [-1, yDoor0, zDoor0], 'x1', [-1, yDoor1, zDoor1], 'roomId1', i, 'roomId2', j, 'wallId1', p, 'wallId2', q);
                        doorList{doorCount} = doorINFO;
                        doorCount = doorCount + 1;   
                        end
                    end               
                    end
                end
                        
            
          
                 % ydir
                % first left second right
                if label0(p) == -2 & label1(q) == 2 & min(wall1(q,1), wall1(q+1,1)) > max(wall0(p,1),wall0(p+1,1)) 
                    % extract wall profile
                  
                    pw1 = floor(mean([x0(p) x0(p+1)]));
                    pw2 = floor(mean([x1(q) x1(q+1)]));
                    w = min(door_param.wpmargin, abs(floor(0.5*(pw1-pw2))));   
                    
                    % first
                    wallProfile0 = zeros(voxel_param.gNum(3), voxel_param.gNum(2));
                    freeSpaceProfile0 = ones(voxel_param.gNum(3), voxel_param.gNum(2));
                    ymin = floor(min(y0(p), y0(p+1)))-door_param.addmargin;
                    ymax = floor(max(y0(p), y0(p+1)))+door_param.addmargin;
                    ymin = max(ymin,1);
                    ymax = min(ymax,voxel_param.gNum(2));                    
                    xx = floor(mean([x0(p) x0(p+1)]));
                    
                    for tt = -door_param.wpmargin : w                                                
                          t = [ymin:ymax];
                          [Z, Y] = meshgrid([1:voxel_param.gNum(3)], t);                 
                          id = (Z-1).*voxel_param.gNum(2)*voxel_param.gNum(1) + (xx+tt-1)*voxel_param.gNum(2) + Y;
                          f2d = zeros(size(freeSpaceProfile0));
                          f2d((Y-1).*voxel_param.gNum(3)+Z) = f3d(id);
                          f2d = imclose(f2d, strel('disk', 5)); 
                          
                          p2d = zeros(size(wallProfile0));
                          p2d((Y-1).*voxel_param.gNum(3)+Z) = p3d(id);                          
                          f2d = imclose(f2d, strel('disk', 5)); 
                          
                          freeSpaceProfile0 = min(freeSpaceProfile0, f2d);
                          wallProfile0 =  wallProfile0 + p2d;
                    end 
                    wallProfile0 = normalize_evidence(wallProfile0, 0);
                    freeSpaceProfile0 = normalize_evidence(freeSpaceProfile0, 0);
                    
                    
                    % second
                    wallProfile1 = zeros(voxel_param.gNum(3), voxel_param.gNum(2));
                    freeSpaceProfile1 = ones(voxel_param.gNum(3), voxel_param.gNum(2));
                    ymin = floor(min(y1(q), y1(q+1)))-door_param.addmargin;
                    ymax = floor(max(y1(q), y1(q+1)))+door_param.addmargin;
                    ymin = max(ymin,1);
                    ymax = min(ymax,voxel_param.gNum(2));                    
                    xx = floor(mean([x1(q) x1(q+1)]));
                    
                    for tt = -w : door_param.wpmargin;                                                
                          t = [ymin:ymax];
                          [Z, Y] = meshgrid([1:voxel_param.gNum(3)], t);                 
                          id = (Z-1).*voxel_param.gNum(2)*voxel_param.gNum(1) + (xx+tt-1)*voxel_param.gNum(2) + Y;
                          f2d = zeros(size(freeSpaceProfile1));
                          f2d((Y-1).*voxel_param.gNum(3)+Z) = f3d(id);

                          p2d = zeros(size(wallProfile1));
                          p2d((Y-1).*voxel_param.gNum(3)+Z) = p3d(id);                          
                          f2d = imclose(f2d, strel('disk', 5));                           
                          freeSpaceProfile1 = min(freeSpaceProfile1, f2d);
                          wallProfile1 =  wallProfile1 + p2d;
                    end 
                    wallProfile1 = normalize_evidence(wallProfile1, 0);
                    freeSpaceProfile1 = normalize_evidence(freeSpaceProfile1, 0);  
                    
                    if abs(pw1-pw2) < door_param.min_walldist/options.leafsize(1)
                    min_x = max(min(y0(p), y0(p+1)),  min(y1(q), y1(q+1)));
                    max_x = min(max(y0(p), y0(p+1)),  max(y1(q), y1(q+1)));
                    [doorY0, doorY1, minCost, flag] = wall_profile_analysis(freeSpaceProfile0, freeSpaceProfile1, wallProfile0, wallProfile1, min_x, max_x, door_param.thresh/options.leafsize(1), i, j, p, q);

                    if flag ~= 0
                        mind = max(abs(doorY0(1) - min_x), abs(doorY1(1) - max_x));
                        figure('WindowStyle', 'docked', 'Name', sprintf('Door is between [Room %d Wall %d] [Room %d Wall %d]', i,p,j,q)),
                        imagesc(point_evidence), axis equal
                        hold on
                        plot([x0(p);x0(p+1)],[y0(p);y0(p+1)],'o-r', 'LineWidth',3)     
                        plot([x1(q);x1(q+1)],[y1(q);y1(q+1)],'o-g', 'LineWidth',3)      
                        title(sprintf('%f %d %d %f %f %f %f', mind, label0(p), label1(q), wall0(p,1), wall0(p,2), wall1(q,1), wall1(q,2)))
                        hold off
                        axis equal
                        
                        yDoor0 = voxel_param.bSize(2)*doorY0(1)/voxel_param.gNum(2)+voxel_param.b0(2);
                        zDoor0 = voxel_param.bSize(3)*doorY0(2)/voxel_param.gNum(3)+voxel_param.b0(3); 
                        yDoor1 = voxel_param.bSize(2)*doorY1(1)/voxel_param.gNum(2)+voxel_param.b0(2);
                        zDoor1 = voxel_param.bSize(3)*doorY1(2)/voxel_param.gNum(3)+voxel_param.b0(3); 
                        xDoor0 = mean([wall0(p,1), wall0(p+1,1)]); 
                        xDoor1 = mean([wall1(q,1), wall1(q+1,1)]);           
                                             
                        aisleFlag = 0;
                        if flag == 2
                           aisleFlag = 1; 
                        end
                        
                        if flag > 0
                        doorINFO = struct('aisleFlag', aisleFlag, 'x0', [-1, yDoor0, zDoor0], 'x1', [-1, yDoor1, zDoor1], 'roomId1', i, 'roomId2', j, 'wallId1', p, 'wallId2', q);
                        doorList{doorCount} = doorINFO;
                        doorCount = doorCount + 1;   
                        end
                    end
                    end           

                end
                
            end
        end

       
        end
        
        
    
         
    end
end
      


doorList = doorList';

end

