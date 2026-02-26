%% params to switch often
loadSkels=1;
debug=1;
lotsaplots=1;
allGraphs=0;

%% model params, taken from previous script
resistivity = 2000;
Rm = 2000;
radius = .01;
Vo = 1;
a=radius/10000;
%x=d/10000;
ra=resistivity/(pi*a^2);
rm=Rm/(2*pi*a);
lc=sqrt(rm/ra);

%% set up locations
baseDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\';
WPN=[baseDir 'Analysis\'];
MPN=[baseDir 'Merge\'];
SMDir=[baseDir 'Analysis\SMs\'];
load([MPN 'obI.mat']);
load([MPN 'dsObj.mat']);
load([WPN 'tis.mat']);

%% load skeletons if necessary
if loadSkels==1
vgcCidList=[2 3 4 5 13];
SMs={};
uniqueVGCs=unique(vgcCidList);
for curSMIt=1:length(uniqueVGCs)
    fileStr=['sm_cid' + string(uniqueVGCs(curSMIt)) + '.mat'];
    if isfile([SMDir + fileStr])
        load([SMDir + fileStr]);
        SMs{curSMIt}=sm;
    end
end
end

%% get the skeleton data into the supermatrix
allEdges=[];
allPos=[];
allSynType=[];
axisObIDList=[];

for curVGC=1:length(SMs)
    curSM=SMs{curVGC};
    distMats{curVGC}=curSM.syn2Skel.syn2SynDist;
    axisObIDList=[axisObIDList;curSM.syn.obID];
    allEdges=[allEdges; curSM.syn.edges];
    allPos=[allPos; curSM.syn.pos];
    allSynType=[allSynType; curSM.syn.synType];
    %allSynType
end

allPreTypes=cid2type(allEdges(:,2),tis);
allPostTypes=cid2type(allEdges(:,1),tis);

bigMat=blkdiag(distMats{1},distMats{2},distMats{3},distMats{4},distMats{5});
infMat=Vo*exp(-bigMat/lc/10000);
infMat(bigMat==0)=0;
typeMat=allSynType*allSynType';
