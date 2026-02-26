function[checkCol  useCol] = getList_seedColor(seedList,checkList)


load('MPN.mat')
load([MPN 'obI.mat']);

seedList = [ 108  201 903 907];

useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
seedPref = seedPreferences(seedList,useList);
allEdges = obI.nameProps.edges(:,[2 1]);




% comb = {1; 2 ; 3 ; 4; ...
%     [1 2]; [1 3] ; [1 4]; [2 3]; ...
%     [2 4]; [3 4]; [1 2 3]; [2 3 4]; ...
%     [1 2 4]; [1 3 4]; [1 2 3 4]}
% 
% col = [1 0 0; 0 1 0;  0 .5 1.5; .5 0 1.5 ; ...
%     1 1 0; 1 .5 1.5; 1 0 1.5;  0 1 1.5; ...
%     .5 1 1.5; .5 .5 1.5; 1 1.5 1.5; .5 1.5 2; ...
%     1.5 1 1.5; 1.5 .5 2; 1 1 1]


comb = {1; 2 ; 3 ; 4; ...
    [1 2]; [1 3] ; [1 4]; [2 3]; ...
    [2 4]; [1 2 3]; [1 2 4]};
     

col = [1 0 0; 0 1 0;  0 .5 1; .5 0 1 ; ...
    1 1 0; 1 .5 1; 1 0 1;  0 1 1; ...
    .2 .9 1 ; 1 1 1; 1 1 1 ];

checkCol = zeros(length(checkList),3);

recC = zeros(length(comb),1);

for i = 1:length(checkList)
    targ = find(seedPref.cellList == checkList(i));
    if isempty(targ)
        seedLink = [];
       for s = 1:length(seedList)
        if sum((allEdges(:,1) == checkList(i)) & (allEdges(:,2) == seedList(s)));
           seedLink = [seedLink s];
        end
       end
        
    else
    seedLink = find(seedPref.sharedSyn(:,targ));
    seedLink = seedLink(:)';
    end
    
    
    
    
    for c = 1:length(comb)
       if length(comb{c}) == length(seedLink);
        if ~sum((comb{c} == seedLink) ==0)
           checkCol(i,:) = col(c,:);
           recC(c) = 1;
           break
        end
       end
    end
    
end

useCol.comb = comb(recC>0);
useCol.col = col(recC>0,:);














