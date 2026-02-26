%function[] = showSM_01(cid)

cid = 0
%makeSM(cid);
load('WPN.mat')

fileName = sprintf('sm_cid%d.mat',cid);
if ~exist([WPN 'SMs\' fileName],'file')
    sm = makeSM(cid);
else
    load([WPN 'SMs\' fileName],'sm')
end

%showRads3D(sm.arbor); 
%showArbor3D(arbor);
%showRadSurf(pos,edge,rad,nodeCol)
%showArborChunks3D(arbor);
%showArborVox(sm.arbor)

if 0
    clf
    renderFVvox(sm.show.voxFV,[0 0 1],.3)
    hold on
    %showRadSurf(sm.nep.pos,sm.nep.edges,sm.nep.nodeRad,[0 1 0],.4)
    showRadSurf(sm.nep.pos,sm.nep.edges,sm.nep.meanNodeRad,[1 0 0],.1)
    showRadSurf(sm.nep.pos,sm.nep.edges,sm.nep.meanNodeRad * 0 + .03,[0 1 0],1)
    hold off
end


if 0  %Show property
    clf
    cSize = 1000;
    cmap = jet(cSize);
    seedPos = sm.nep.pos(sm.nep.seedNode,:);
    col = driftNetCol(sm.nep,.02);
%     dists = getDist(sm.nep.pos,seedPos);
%     colInd = floor(dists/max(dists) * (cSize-1)) + 1;
    showRadSurf(sm.nep.pos,sm.nep.edges,sm.nep.meanNodeRad * 0 + .03,col,1);
    scatter3(seedPos(1),seedPos(2),seedPos(3),300,'o','filled','w')
    hold off
end


if 1
    clf
    renderFVvox(sm.show.voxFV,[0 0 1],.4)
    hold on
    %showRadSurf(sm.nep.pos,sm.nep.edges,sm.nep.nodeRad,[0 1 0],.4)
%   showRadSurf(sm.nep.pos,sm.nep.edges,sm.nep.meanNodeRad,[1 0 0],.1)
    showRadSurf(sm.nep.pos,sm.nep.edges,sm.nep.meanNodeRad * 0 + .03,[1 1 0],1)
    hold off
end



%% Render sub
if 0
    downSamp = 1;
    renderProps.smooth = 0;
    renderProps.resize = 2;
    renderProps.smoothPatch = 10;
    fv = subVolFV(smallSub,[],renderProps);
    [p] = renderFV(fv,col(i,:),alph(i));
    view([0 0])
    axis off
    pause(.01)
    hold on
end


%% Test smooth patch
if 0
    clf
    fv = showRadSurf(sm.nep.pos,sm.nep.edges,sm.nep.meanNodeRad,[1 0 0],.2)
    clf
    fv3 = face4to3(fv);
    sMode = 1;
    lambda = .8;
    sigma = 1;
    itt = 2;
    fv3s=smoothpatch(fv3,sMode,itt,lambda,sigma);
    clf
    %[p] = renderFV(fv,[1 0 0],.2);
    [p] = renderFV(fv3,[0 1 0],.2);
    [p] = renderFV(fv3s,[0 0 1],.5);
end



            
            
            
            
            
            
            
            
            
            
            
            
            



