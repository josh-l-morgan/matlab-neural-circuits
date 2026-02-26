function output=getDepthDat(z,x,y,GCL,INL,histBinEdges)
    [zg,zl,iplDepths]=getIPLdepth(z,x,y,GCL,INL);
    depthHistDat=histcounts(iplDepths,'BinEdges',histBinEdges);
    %get a common space of 0-1 with 0.01 bins for plotting all data.
    arborDepthDat=zeros(length(depthHistDat),2);
    arborDepthDat(:,2)=depthHistDat;
    arborDepthDat(:,1)=histBinEdges(1:end-1)+(histBinEdges(2)-histBinEdges(1))/2;
    output=arborDepthDat;
end