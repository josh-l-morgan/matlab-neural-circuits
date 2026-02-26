function[swc] = swc2fv(SPN)

global glob globSWC

dSPN = dir([SPN '*.swc']);
sNams = {dSPN.name};

libDir = globSWC.fvDir;

%% make swc
cellNum = length(sNams);
swc.cells.cids = 1:cellNum;
swc.cells.label = sNams;
swc.cells.type.typeID = zeros(cellNum,1); 

%% run each file in the folder
tempFig = figure;
allSyn.pos = [];
allSyn.pre = [];
allSyn.post = [];
for i = 1:length(sNams)
    nam = sNams{i};
     cid = i;
     
    fileName = [SPN nam];
    
    [nep, syn] = swc2nep(fileName);
    
    allSyn.pos = cat(1,allSyn.pos,syn.pos);
    allSyn.pre = cat(1,allSyn.pre,syn.pre * cid);
    allSyn.post = cat(1,allSyn.post,syn.post * cid);
    
    dif = nep.pos(nep.edges(:,1),:) - nep.pos(nep.edges(:,2),:);
    dist = sqrt(sum(dif.^2,2));
    minNodeRad = median(dist)/10;
    
    
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
   
    fvFilename = sprintf('%s%d.mat',libDir,cid);
    save(fvFilename,'fv')
    
    
end
close(tempFig)


%Combine syns

pos = allSyn.pos;
maxPos = max(pos,[],1);
ind = sub2ind(maxPos,pos(:,1),pos(:,2),pos(:,3));
[uInd ia ib] = unique(ind);

isPre = find(syn.post == 0);
isPost = find(syn.pre == 0);
openPost = isPost*0+1;

dif1 = pos(isPre,1) - pos(isPost,1)';
dif2 = pos(isPre,2) - pos(isPost,2)';
dif3 = pos(isPre,3) - pos(isPost,3)';
dist = sqrt(dif1.^2 + dif2.^2 + dif3.^2);



for i = 1:length(isPre)
    
    
    
end




swc.syn = allSyn;
