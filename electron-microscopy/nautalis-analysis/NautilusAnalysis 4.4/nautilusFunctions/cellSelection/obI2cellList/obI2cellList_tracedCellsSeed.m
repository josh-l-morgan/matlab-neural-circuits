function[useList] = obI2cellList_tracedCells(obI,seedList);

%%Make list of all traced thalamocortical cells and the rgcs that innervate
%%them

tracedTCR = [   108   109   116   117   120   123   129   133   134   148   156   159   162 ... 
   163   169   170   201   203   205   206   207   210   212   213   215   216 ...
   217   218   232   267   601   907   110   268   273   237   224];


%% subjectively chosen to have substantial tracing performed
 tracedAxons = [ 1001        1006        1009        1012        1014        1021        1023        1025        1027        1028        1029        1030        1031 ...
        1032        1033        1034        1036        1037        1040        1041        1050        1051        1053        1054        1055        1056 ...
        1058        1059        1060        2003        2004        2006        2007        2008        2009        2011        2016        2024        2027 ...
        2028        2030        2032        2033        2034        2035        2038        2041        5003        5004        5005        5102        5104 ...
        5106        5108        6000        6001        6101        7001        7004        7016        7018        7019        7021        7022        9101  9102];

%%


if ~exist('seedList','var')
    seedList = [];
end

rgcList = obI.nameProps.cellNum(obI.nameProps.rgc);
rgcList = unique(rgcList(rgcList>0));


tcrList = obI.nameProps.cellNum(obI.nameProps.tcr);
tcrList = unique(tcrList(tcrList>0));


%cellList =  tracedTCR;

synapses = obI.nameProps.edges;
edges = synapses(:,1:2);



isPre = [];
for i = 1:length(seedList)
    targ = seedList(i);
    isTarg = synapses(:,1) == targ;
    isPre = [isPre synapses(isTarg,2)'];
end

%}
axList = intersect(tracedAxons,rgcList);
axList = intersect(axList,isPre);

isPost = [];
for i = 1:length(axList)
    targ = axList(i);
    isTarg = synapses(:,2) == targ;
    isPost = [isPost synapses(isTarg,1)'];
end

cellList = intersect(tracedTCR,tcrList);
cellList = intersect(cellList,isPost);





%% graph

con = zeros(length(axList),length(cellList));
for i = 1:length(axList)
    for p = 1:length(cellList)
        con(i,p) = sum( (edges(:,1) == cellList(p)) & (edges(:,2) == axList(i)));
    end
end


%% Create outpute structure

useList.seedList = seedList;
useList.preList = axList;
useList.postList = cellList;
useList.con = con;




