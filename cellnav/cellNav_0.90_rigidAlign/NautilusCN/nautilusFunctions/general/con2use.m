function[useList] = con2use(dat)

 useList.con = dat ;
    useList.preList = 1:size(dat,1);
    useList.postList = size(dat,1)+1: size(dat,1)+size(dat,2);
    useList.seedList = [];
    seedList = [];
    
    allEdges = con2syn(dat);
    
    allEdges(:,1) = useList.preList(allEdges(:,1));
    allEdges(:,2) = useList.postList(allEdges(:,2));
    useList.allEdges = allEdges;
    
    
    [y x] = find(dat);
    allWeights = [y x dat(dat>0)];
    useList.allWeights = allWeights;
    