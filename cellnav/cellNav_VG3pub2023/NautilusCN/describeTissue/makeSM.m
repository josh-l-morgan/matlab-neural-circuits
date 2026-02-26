function[sm] = makeSM(targID)


load('MPN.mat')
load([MPN 'obI.mat'])
load([MPN 'dsObj.mat'])
%load('WPN.mat');
load([WPN 'tis.mat']);

%% Cell metadata
sm.cid = targID;
sm.fileName = sprintf('sm_cid%d.mat',sm.cid);
targ = find(tis.cids==sm.cid);
sm.cell.anchor = tis.cells.anchors(targ,:);
sm.cell.typeID = tis.cells.type.typeID(targ);
sm.cell.subTypeID = tis.cells.type.subTypeID(targ);

%% Synapses
sm = getSynSM(sm);

%% Skeleton
%sm = makeSkeletonForSM(sm); % skeletonize the arbor
sm = makeBranchedSkeletonForSM(sm); % skeletonize the arbor
%showRadSurf(sm.nep.pos,sm.nep.edges,sm.nep.nodeRad)
if ~exist([WPN 'SMs\'],'dir'),mkdir([WPN 'SMs\']);end 
save([WPN 'SMs\' sm.fileName],'sm','-v7.3')

if ~isempty(sm.nep)
    %sm = getSkelForSM(sm,.5); % up sample edges so that there can be multiple nodes for each current edge
    sm = skel2skelDist(sm);
    save([WPN 'SMs\' sm.fileName],'sm','-v7.3')
    sm = getSkelProps(sm);
    
    sm.show.voxFV = sm.arbor.vox.fv;
    vert = sm.show.voxFV.vertices;
    vert = vert .* repmat(sm.arbor.voxSize,[size(vert,1) 1]);
    sm.show.voxFV.vertices = vert;
    % sm.show.voxFV.vertices = vert+ repmat(sm.arbor.offset .* sm.arbor.voxSize,...
    

    
    sm = syn2SkelSM(sm);
    %save([WPN 'SMs\' sm.fileName],'sm','-v7.3')      
  
   
    sm = sm2swc(sm) %write swc file
    sm.nep.swcS = nep2swc(sm.nep);
    save([WPN 'SMs\' sm.fileName],'sm','-v7.3')
    
    
    
    
    %% Show result
    clf
    renderFVvox(sm.show.voxFV,[0 0 1],.4)
    hold on
    showRadSurf(sm.nep.pos,sm.nep.edges,sm.nep.meanNodeRad,[0 1 0],.4)
    %   showRadSurf(sm.nep.pos,sm.nep.edges,sm.nep.meanNodeRad,[1 0 0],.1)
    %   showRadSurf(sm.nep.pos,sm.nep.edges,sm.nep.meanNodeRad * 0 + .03,[1 1 0],1)
    hold off
    
end


% %% Save Nep
% nepName = sprintf('nep_cid%d.mat',sm.cid);
% nep = sm.nep;
% if ~exist([WPN 'neps\'],'dir'),mkdir([WPN 'neps\']),end
% save([WPN 'neps\' nepName],'nep','-v7.3')

smxName = sprintf('smx_cid%d.mat',sm.cid);
save([WPN 'SMs\' smxName],'-struct','sm','nep','syn',...
    'cell','skel2skel','syn2Skel','-v7.3');

%% Old

%sm = addDatToSynMat(obI,targID);
%sm = getTopoEucDistBetweenSyn(sm);
%sm = labelShaftSkel(sm);
%sm = labelSubTypes(sm);
%sm = synSkel2synSkelDist(sm);











