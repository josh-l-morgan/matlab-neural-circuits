%synDoubleChecker

load('Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Analysis\tis.mat');

synEucDistMat=zeros(length(tis.syn.pos),length(tis.syn.pos));
for synIt=1:length(tis.syn.pos)
    curSynPos=tis.syn.pos(synIt,:);
    curSynDisps=tis.syn.pos-curSynPos;
    curSynDists=sqrt(sum(curSynDisps.*curSynDisps,2));
    synEucDistMat(:,synIt)=curSynDists;
end

%get the synapses with no anchor point
cardinalSyns=tis.syn.obID(tis.syn.pos(:,1)<1);
badSynIDs=find(tis.syn.pos(:,1)<1);
badNames=tis.obI.nameProps.names(cardinalSyns);


