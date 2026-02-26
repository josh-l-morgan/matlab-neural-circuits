function[] = swc2skel(SPN)

global glob

dSPN = dir([SPN '*.swc']);
sNams = {dSPN.name};

libDir = glob.fvDir;


%% run each file in the folder
for i = 1:length(sNams)
    nam = sNams{i};
    fileName = [SPN nam];
    
    nep = swc2nep(fileName)
    
    
    dif = nep.pos(nep.edges(:,1),:) - nep.pos(nep.edges(:,2),:);
    dist = sqrt(sum(dif.^2,2));
    minNodeRad = median(dist)/10;
    
    tempFig = figure;
    nodeRad = nep.nodeRad;
    %nodeRad(nodeRad<minNodeRad) = minNodeRad;
    nodeRad(:) = minNodeRad;
    fv = showRadSurf_cnv(nep.pos,nep.edges,nodeRad,[1 1 1],1,tempFig);
    
    
    
    %% nodes
    isRad = find(nep.nodeRad>0);
    [x,y,z] = sphere;
    fvN = surf2patch(x,y,z);
    svFac = fvN.faces;
    svVert = fvN.vertices;
    
    for n = 1:length(isRad)
        %surf(x,y,z)
        numVert = size(fv.vertices,1);
        
        addVert = svVert * nep.nodeRad(isRad(n)) + repmat(nep.pos(n,:),[size(svVert,1) 1]);
        addFac = svFac + numVert;
        addCol = ones(size(svVert,1),3);
        
        fv.vertices = cat(1,fv.vertices,addVert);
        fv.faces = cat(1,fv.faces,addFac);
        fv.FaceVertexCData = cat(1,fv.FaceVertexCData,addCol);
        
        
    end
    patch(fv)
    pause(.01)
    cid = i;
    fvFilename = sprintf('%s%d.mat',libDir,cid);
    save(fvFilename,'fv')
    close(tempFig)
    
    
end


