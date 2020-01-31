function convert_floorplan_to_ply(Path)
% READ_FLOORPLAN - convert a floorplan.txt to ply
% READ_FLOORPLAN(Path): 
%                          Path: the filename of the floorplan.txt
%                          OUTPUT: a ply file named '[Path].ply'

[fid,Msg] = fopen(Path,'rt');	% open file in read text mode
if fid == -1, error(Msg); end

% read transformation matrix
rotation = (fscanf(fid, '%f', [3 3]))';
numRooms = fscanf(fid, '%d', 1);

withFloor = 1;
withCeil = 0;

PTS_all = cell(numRooms,1);
TRI_all = cell(numRooms,1);
numPtsRoomWall = cell(numRooms,1);

numPtsTotal = 0;
WALLS = cell(numRooms, 1);
for i = 1 : numRooms
    numTags = fscanf(fid, '%d', 1);
    tags = cell(numTags, 1);
    
    % read tags
    for j = 1 : numTags
        tags{j,1} = fscanf(fid, '%s', 1);
    end
    
    numWalls = fscanf(fid, '%d', 1);
    walls = (fscanf(fid, '%f', [2 numWalls]))';
    walls(end+1,:) = walls(1,:);
    WALLS{i, 1} = walls;
    fH = fscanf(fid, '%f', 1);
    cH = fscanf(fid, '%f', 1);
    
    % read walls
    PTS_room = cell(numWalls+2,1);
    TRI_room = cell(numWalls+2,1);
    numPts = 0;
    numPtsRoom = zeros(numWalls,1);
    for j = 1 : numWalls
        numVertices = fscanf(fid, '%d', 1);
        numTrianglesWall = fscanf(fid, '%d', 1);

        pts = (fscanf(fid, '%f', [2 numVertices]))';
        tri = (fscanf(fid, '%f', [10 numTrianglesWall]))';
        X = zeros(numVertices, 3);

        for k = 1 : numVertices
           X(k, :) = [walls(j,:) + pts(k,1)*(walls(j+1,:)-walls(j,:)), fH + pts(k,2)*(cH-fH)];              
        end

        PTS_room{j,1} = X;
        numPtsRoom(j) = numPtsTotal + numPts;
        TRI_room{j,1} = tri(:,1:3) + numPts;
        numPts = numPts + numVertices;
    end     
    numPtsRoomWall{i} = numPtsRoom;
    
    % read floors and ceils    
    Xfloor = zeros(numWalls, 3);
    Xceil = zeros(numWalls, 3);
    for j = 1 : numWalls
         Xfloor(j,:) = [walls(j,:), fH]; % wall-floor corners
         Xceil(j,:) = [walls(j,:), cH]; % wall-floor corners
    end   
    
    numTrianglesFloor = fscanf(fid, '%d', 1);    
    tri = (fscanf(fid, '%f', [10 numTrianglesFloor]))';
    count = 1;
    if withFloor
    PTS_room{numWalls+count, :} = Xfloor;     
    TRI_room{numWalls+count,1} = tri(:, 1:3) + numPts;
    numPts = numPts + numWalls;
    count = count + 1;
    end
    
    numTrianglesCeil = fscanf(fid, '%d', 1);
    tri = (fscanf(fid, '%f', [10 numTrianglesCeil]))';
    if withCeil    
    PTS_room{numWalls+count, :} = Xceil;
    TRI_room{numWalls+count,1} = tri(:, 1:3) + numPts;
    numPts = numPts + numWalls;   
    end
    
    PTS_all{i, 1} = cell2mat(PTS_room);
    TRI_all{i, 1} = cell2mat(TRI_room) + numPtsTotal;
    numPtsTotal = numPtsTotal + numPts;
    

end


% read door connection 
numDoorConnection = fscanf(fid, '%d', 1);

TRI_door = cell(numDoorConnection,1);

for i = 1 : numDoorConnection
    
    roomid1 = fscanf(fid, '%d', 1); % 1st room id that the door is in
    wallid1 = fscanf(fid, '%d', 1); % 1st wall id that the door is on
    roomid2 = fscanf(fid, '%d', 1); % 2nd room id that the door is in
    wallid2 = fscanf(fid, '%d', 1); % 2nd wall id that the door is on
    didx1 = fscanf(fid, '%d', 4); % point index (1st wall)
    didx2 = fscanf(fid, '%d', 4); % point index (2nd wall)
    triconnection = (fscanf(fid, '%d', [10 8]))';
    triconnection = triconnection(:,1:3);
    didx1 = didx1 + numPtsRoomWall{roomid1+1}(wallid1+1); % 0 1 2 3
    didx2 = didx2 + numPtsRoomWall{roomid2+1}(wallid2+1); % 4 5 6 7
    
    source_idx = [0:7]';
    target_idx = [didx1(:);didx2(:)];
    triconnection_transformed = zeros(size(triconnection));
    for k = 1:8
        triconnection_transformed(find(triconnection==source_idx(k))) = target_idx(k);
    end
    TRI_door{i} = triconnection_transformed;

end



figure, 
TRI = [cell2mat(TRI_all);cell2mat(TRI_door)];
PTS = cell2mat(PTS_all);
trimesh(TRI+1,PTS(:,1),PTS(:,2),PTS(:,3)) 
title('Recovered geometry w/o ceiling');
axis equal



plyname = sprintf('%s.ply', Path);
if ~isempty(plyname)
simple_write_mesh_ply((rotation*PTS')', TRI, plyname);
end




