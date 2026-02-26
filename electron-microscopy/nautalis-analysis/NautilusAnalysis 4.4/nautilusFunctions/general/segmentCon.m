function[cellGroup, groupNum] = segmentCon(con);

if 0 %Detect
    %% detect directed
    directed = 1;
    [y x] = size(con);
    if y == x
        if  ~sum(con-con');
            directed = 0;
        end
    end
    
    
    if directed  %% convert to undirected symmetric graph
        
        con2 = zeros(y+x);
        con2(1:y,y+1:end) = con>0;
        con2(y+1:end,1:y) = con'>0;
        %image(con2)
    else
        
        con2 = con;
    end
else
    con2 = con + con'; %force symmetric
end


num = size(con2,1); %node number
cellGroup = zeros(1,num); %list of groups for nodes
for g = 1:num; % run number of groups less then or equal to node number
    
    
    nextY = find(cellGroup ==0,1); % first node with no group
    while 1 % always
        Y = unique(nextY); %start using found nodes
        cellGroup(Y) = g; %assign found nodes to group number
        con2(:,Y) = 0; %remove current nodes for receiving connections
        [oldY nextY] = find(con2(Y,:)>0); 
        if isempty(nextY),break,end
    end
    if ~sum(cellGroup==0),break,end
end

groupNum = g;
