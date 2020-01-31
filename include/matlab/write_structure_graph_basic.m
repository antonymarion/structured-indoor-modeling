
function write_structure_graph_basic(filename, WALLS, LABELS, doorList, floorHeight, ceilHeight, rotmat)
% WRITE_STRUCTURE_GRAPH_BASIC(filename, WALLS, LABELS, doorList, floorHeight, ceilHeight, rotmat)
% Export the structure graph about the basic architectural elements ( room,
% wall, floor, ceiling, door)
% NOTE the wall direction is clockwise 
% NOTE the polygon direction is counter-clockwise (frontal face of a wall is pointing to camera)

display('Export structure graph in .txt format...');
[fid, Msg] = fopen(filename, 'wt');
fprintf(fid, '%g %g %g\n', rotmat);

numRooms = length(WALLS); 

if isempty(cell2mat(doorList));
    numDoors = 0;
else
    numDoors = length(doorList);
end
fprintf(fid, '%d\n', numRooms); % number of rooms

doorIndices1 = cell(length(doorList),1);
doorIndices2 = cell(length(doorList),1);

%% for each room
 numCurTri = 0;
for i = 1 : length(WALLS)
    words = {['bedroom'];['showerroom']};
    numWords = length(words);
    walls = WALLS{i};
    labels = LABELS{i};
    numWalls = size(WALLS{i},1);    
    
    fprintf(fid, '%d\n', numWords); % number of words describing the room type
    for k = 1 : numWords
    fprintf(fid, '%s\n', words{k}); % words
    end     
    fprintf(fid, '%d\n', numWalls); % number of the following 2d coordinatges
    fprintf(fid, '%g %g\n', walls'); % 2d coordinates in the manhattan coordinate frame
    fprintf(fid, '%d\n', floorHeight(i)); % floor_height
    fprintf(fid, '%d\n', ceilHeight(i)); % ceiling_height
    
    doorWallId = [];
    doorRoomId = [];
    doorID = [];
   
    for k = 1 : numDoors
        if doorList{k}.roomId1 == i
            doorWallId(end+1) = doorList{k}.wallId1;  
            doorID(end+1) = k;
            doorRoomId(end+1) = 1;
        end

        if doorList{k}.roomId2 == i
            doorWallId(end+1) = doorList{k}.wallId2;   
            doorID(end+1) = k;
            doorRoomId(end+1) = 2;
        end
    end
    walls(end+1,:) = walls(1,:);
    % for each wall  

    for k = 1 : numWalls
        uvlist = [];
        trilist = [];
        doorOnThisWall = unique(find(doorWallId==k));
        if isempty(doorOnThisWall) % walls w/o door
            uv = [0.0 0.0;1.0 0.0;1.0 1.0;0.0 1.0];
            uvlist(end+1:end+size(uv,1), :) = uv;
            tri = [0 2 3;0 1 2];
            if ~isempty(trilist)
            trilist(end+1:end+size(tri,1), :) = tri;  
            else
            trilist = tri;  
            end
                
        else % walls w/ door
           uv = [0.0 0.0;1.0 0.0;1.0 1.0;0.0 1.0];
           uvlist(end+1:end+size(uv,1), :) = uv;
           
           
           uvlist_temp = cell(length(doorOnThisWall),1);
            %% add new verteces
            
        
           w0 = walls(k,:);
           w1 = walls(k+1,:);
           dirwall = sign(sum(w1-w0));
%            dirwall = sign(w1(1)-w0(1));
           
           for j = 1 : length(doorOnThisWall)
               didx = doorID(doorOnThisWall(j));
               ridx = doorRoomId(doorOnThisWall(j));
               
               
               if labels(k) == -1 | labels(k) == 1 % wall is horizontal to the MW coordinate
                   x0 = [doorList{didx}.x0(1), doorList{didx}.x0(3)];
                   x1 = [doorList{didx}.x1(1), doorList{didx}.x1(3)];                   
                   dirdoor = sign(x1(1)-x0(1));
                   if dirwall ~= dirdoor
                       x0_temp = x0;
                       x0 = x1;
                       x1 = x0_temp;
                       a = x0(2);
                       x0(2) = x1(2);
                       x1(2) = a;
                   end


                   u0 = (x0(1) - w0(1))/(w1(1)-w0(1));
                   v0 = (x0(2) - floorHeight(i))/(ceilHeight(i)-floorHeight(i));
                   u1 = (x1(1) - w0(1))/(w1(1)-w0(1));
                   v1 = (x1(2) - floorHeight(i))/(ceilHeight(i)-floorHeight(i));
                   uv = [u0 0.0;u1 0.0;u1 v1;u1 1.0;u0 1.0;u0 v1];                                    
                   uvlist_temp{j,1} = uv;
                   
               
               else
                   x0 = [doorList{didx}.x0(2), doorList{didx}.x0(3)];
                   x1 = [doorList{didx}.x1(2), doorList{didx}.x1(3)];                                                       
                   dirwall = sign(w1(2)-w0(2));
                   dirdoor = sign(x1(1)-x0(1));
                   
                   if dirwall ~= dirdoor
                       x0_temp = x0;
                       x0 = x1;
                       x1 = x0_temp;
                       a = x0(2);
                       x0(2) = x1(2);
                       x1(2) = a;
                   end
              
                   u0 = (x0(1) - w0(2))/(w1(2)-w0(2));
                   v0 = (x0(2) - floorHeight(i))/(ceilHeight(i)-floorHeight(i));
                   u1 = (x1(1) - w0(2))/(w1(2)-w0(2));
                   v1 = (x1(2) - floorHeight(i))/(ceilHeight(i)-floorHeight(i));
                   uv = [u0 0.0;u1 0.0;u1 v1;u1 1.0;u0 1.0;u0 v1];                  
                   uvlist_temp{j,1} = uv;
                   
               end
           end    

         
           idx_order = ones(length(doorOnThisWall),2);
           if length(uvlist_temp) > 1
               
               uvlist_temp_copy = uvlist_temp;
               order_list = zeros(length(uvlist_temp),1);
               for j = 1 : length(uvlist_temp) 
                   order_list(j) = uvlist_temp{j}(1,1);
               end
               
              
               idx_order = sortrows([[1:length(doorOnThisWall)]',order_list], 2);    

               for j = 1 : length(doorOnThisWall) 
                  uvlist_temp{j} = uvlist_temp_copy{idx_order(j,1),1};
               end
              
           end
           
%            %%%%%%%%%%%%%%%%%%%%%%
%            if length(uvlist_temp) > 1
%                cell2mat(uvlist_temp)
% %                dirwall
% %                idx_order
% %                idx_order = [2 3 1]'
%            end
%            %%%%%%%%%%%%%%%%%%%%%%%%
           
           
           uvlist_temp = cell2mat(uvlist_temp);
           uvlist(end+1:end+size(uvlist_temp,1), :) = uvlist_temp;
           
         %% add new triangles
           for j = 1 : length(doorOnThisWall)
               didx = doorID(doorOnThisWall(j));
               ridx = doorRoomId(doorOnThisWall(j));
               if ridx == 1
                  doorIndices1{didx} = 4 + 6*(find(idx_order(:,1)==j)-1) + [0 1 2 5];
                  
               else
                  doorIndices2{didx} = 4 + 6*(find(idx_order(:,1)==j)-1) + [0 1 2 5]; 
               end
           end
           
           tri = [0 8 3;0 9 8;0 4 9]; % for first door
           for j = 1 : length(doorOnThisWall)-1
               rootIndex = 4 + 6*(j-1); % the left bottom index of j-th door
               tridoor = [5 3 4;5 2 3]; % for triangle above the door
               tri(end+1:end+size(tridoor,1),:) = rootIndex + tridoor;
               tridoor = [2 10 3;1 10 2;1 11 10;1 6 11]; % for triangle between the door
               tri(end+1:end+size(tridoor,1),:) = rootIndex + tridoor;
           end
           
           tri(end+1:end+5,:) = [4 + 6*(length(doorOnThisWall)-1) + 5, 4 + 6*(length(doorOnThisWall)-1) + 3, 4 + 6*(length(doorOnThisWall)-1) + 4;...
               4 + 6*(length(doorOnThisWall)-1) + 5, 4 + 6*(length(doorOnThisWall)-1) + 2, 4 + 6*(length(doorOnThisWall)-1) + 3;...
               4 + 6*(length(doorOnThisWall)-1) + 2, 2, 4 + 6*(length(doorOnThisWall)-1) + 3;...
               4 + 6*(length(doorOnThisWall)-1) + 1, 2, 4 + 6*(length(doorOnThisWall)-1) + 2;...
               4 + 6*(length(doorOnThisWall)-1)+1 1, 2]; % for last door     
                     
            if ~isempty(trilist)
            trilist(end+1:end+size(tri,1), :) = tri;  
            else
            trilist = tri;  
            end            
        end
        
        %% write the wall
        numVertices = size(uvlist,1);
        numTriangle = size(trilist,1);
        fprintf(fid, '%d %d\n', numVertices, numTriangle); % ceiling_height
        fprintf(fid, '%g %g\n', uvlist');
        for j = 1 : size(trilist,1)
        fprintf(fid, '%d %d %d -1 -1 -1 -1 -1 -1 -1\n', trilist(j,1),trilist(j,2),trilist(j,3));
        end
        numCurTri = numCurTri + size(trilist,1);
       

      
        
    end
    %% for floor
    C = zeros(numWalls,2);
    for j = 1 : numWalls - 1
    C(j,1) = j;
    C(j,2) = j+1;
    end
    C(numWalls,1) = numWalls;
    C(numWalls,2) = 1;

    DT = delaunayTriangulation(WALLS{i}, C);
    IO = isInterior(DT);
    numTri = length(find(IO==1));
    triFloor = DT(IO, :)-1;
    fprintf(fid, '%d\n', size(triFloor,1));
    for j = 1 : size(triFloor,1)
        fprintf(fid, '%d %d %d -1 -1 -1 -1 -1 -1 -1\n',triFloor(j,1), triFloor(j,2), triFloor(j,3));
    end   
    fprintf(fid, '%d\n', size(triFloor,1));
    for j = 1 : size(triFloor,1)
        fprintf(fid, '%d %d %d -1 -1 -1 -1 -1 -1 -1\n',triFloor(j,1), triFloor(j,3), triFloor(j,2));
    end    
    
end
    % for door information

    fprintf(fid, '%d\n', numDoors); 
    for j = 1 : numDoors
    didx1 = doorIndices1{j};
    didx2 = doorIndices2{j};
    fprintf(fid, '%d %d %d %d\n', doorList{j}.roomId1-1, doorList{j}.wallId1-1, doorList{j}.roomId2-1, doorList{j}.wallId2-1);
    fprintf(fid, '%d %d %d %d\n', didx1');
    fprintf(fid, '%d %d %d %d\n', didx2');
       
    triConnection = zeros(8, 3);   
    triConnection(1,:) = [0 6 3];   
    triConnection(2,:) = [0 5 6];   
    
    triConnection(3,:) = [4 2 7];   
    triConnection(4,:) = [4 1 2]; 

    triConnection(5,:) = [5 1 4];   
    triConnection(6,:) = [5 0 1];   
    
    triConnection(7,:) = [3 7 2];   
    triConnection(8,:) = [3 6 7];
    for p = 1 : 8
    fprintf(fid, '%d %d %d -1 -1 -1 -1 -1 -1 -1\n', triConnection(p,1), triConnection(p,2), triConnection(p,3));
    end
    
    end

fclose(fid);
    
end

