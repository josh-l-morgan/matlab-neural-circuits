
load('MPN.mat')
load([MPN 'obI.mat']);
seedList = [ 108 201 109 907 903];

[axWithPrim axWithPrim201] = getList_primaries;
tempList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
tracedList = obI2cellList_tracedCells(obI);
glomList = getList_glomID([1,2]);
cellList = tempList.postList;

useList = obI2cellList_prePostList2Use(obI,axWithPrim201,cellList);

seedPref = seedPreferences(seedList,useList);
allEdges = obI.nameProps.edges(:,[2 1]);
con = useList.con;

image(con*20)