
if 0
    clear all
    [TFN TPN] = uigetfile
    load([TPN TFN])
else
    clear all
    TFN = 'sm_cid4.mat';
    TPN = 'G:\IxQ\Matlab\Analysis\SMs\';
    load([TPN TFN])
end



pos = sm.nep.pos;
edges = sm.nep.edges;
groupEdges = sm.nep.groupEdges;
bones = sm.nep.bones;
rad = sm.nep.nodeRad;
rad = sm.nep.meanNodeRad;

sm.nep.props