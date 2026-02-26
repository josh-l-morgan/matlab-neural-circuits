function[useList] = filtUseList(useList,useNodes)

[seedList use] = intersect(useList.seedList, useNodes);
useList.seedList = seedList;
useList.preList = useList.preList(use);
useList.postList = useList.postList(use);
useList.nodes = useList.nodes(use);
useList.nodeType = useList.nodeType(use);
useList.con = useList.con(use,use);