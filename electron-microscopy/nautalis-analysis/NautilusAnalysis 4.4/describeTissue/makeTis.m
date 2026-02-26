%function[tis]  = makeTis


%% load
load('MPN.mat')
%load('WPN.mat')

load([MPN 'obI.mat'])
load([MPN 'dsObj.mat'])
if exist([MPN 'dat.mat'],'file'),load([MPN 'dat.mat'])
dat = parseGoogleDat;
else
    dat.cid = [];
end

%% copy from obI
tis.obI = obI;
tis.cids = obI.cell.name;
tis.cells.cids = obI.cell.name;
tis.cells.mainID = obI.cell.mainID;
tis.cells.label = obI.cell.label;
tis.dat = dat;


%% Cell types
type.typeNames = {'rgc' 'tcr' 'lin' 'unk','phot','hrz','bpc','amc','mul'};
type.typeNameIDs = 1:length(type.typeNames);
numSegs = length(obI.nameProps.names);
typeMat = zeros(numSegs,length(type.typeNames));
for i = 1:length(type.typeNames)
    if isfield(obI.nameProps.tag,type.typeNames{i});
        val = getfield(obI.nameProps.tag,type.typeNames{i});
        typeMat(:,i) = val;
    end
end
type.cellTypeMat = typeMat(obI.cell.mainID,:);

for i = 1:length(tis.cells.cids)
    foundType = find(type.cellTypeMat(i,:));
    
    type.typeID(i) = 0; %by default
    if ~isempty(foundType)
        tis.cells.info(i).typeTags = type.typeNames{foundType};   
    else
        tis.cells.info(i).typeTags = 'none';
    end
    if length(foundType) == 1
        type.typeID(i) = foundType;
    else
        type.typeID(i) = 0;
    end
        
end

for i = 1:length(type.typeNames)
    tag = type.typeNames{i};
    eval(sprintf('typeLists.%s = tis.cells.cids(find(type.cellTypeMat(:,i)));',tag));
end


%% cell Subtypes

type.subTypeNames{8} = {'vgc'};
type.subTypeNames{7} = {'bc1' 'bc2' 'bc3a' 'bc3b' 'bc4' 'bc5i' 'bc5o' 'bc5t' 'bc6' 'bc7' 'bc8' 'bc11' 'bcunk'};


for t = 1:length(type.subTypeNames)
    checkNames = type.subTypeNames{t};
    for s = 1:length(checkNames)
        nam = checkNames{s};
        if isfield(obI.nameProps.tag,checkNames{s});
            val = getfield(obI.nameProps.tag,checkNames{s});
            subTypeMat(t,s,:) = val;
        end
    end
end



for i = 1:length(tis.cells.cids)
   
    subTypeTags = subTypeMat(:,:,tis.cells.mainID(i));
    tis.cells.info(i).subTypeTags = subTypeTags;
    subType = find(sum(subTypeTags,1),1);
    
    [T S] = find(subTypeTags);
    
    
    if isempty(subType)
        type.subTypeID(i) = 0;
    elseif length(T) > 1
        T = T(1);
        S = S(1);
        disp(sprintf('cid %d has multiple subtypes',tis.cells.cids(i)))
    else
        type.typeID(i) = T;
        type.subTypeID(i) = S;
    end
    
end


%{
type.subTypeNameIDs = 1;
numSegs = length(obI.nameProps.names);
subTypeMat = zeros(numSegs,length(type.subTypeNames));
for i = 1:length(type.subTypeNames)
    if isfield(obI.nameProps.tag,type.subTypeNames{i});
        val = getfield(obI.nameProps.tag,type.subTypeNames{i});
        subTypeMat(:,i) = val;
    end
end
type.cellSubTypeMat = subTypeMat(obI.cell.mainID,:);


for i = 1:length(tis.cells.cids)
    foundType = find(type.cellSubTypeMat(i,:));
    
    type.subTypeID(i) = 0; %by default
    if ~isempty(foundType)
        tis.cells.info(i).subTypeTags = type.subTypeNames{foundType};   
    else
        tis.cells.info(i).subTypeTags = 'none';
    end
    if length(foundType) == 1
        type.subTypeID(i) = foundType;
    else
        type.subTypeID(i) = 0;
    end
        
end

for i = 1:length(type.subTypeNames)
    tag = type.subTypeNames{i};
    eval(sprintf('subTypeLists.%s = find(type.cellSubTypeMat(:,i));',tag));
end
%}

%% on off
type.onOff = zeros(length(type.typeID),2);

%% tracing
tis.cells.tracing = zeros(length(tis.cells.cids),1);

%% collect
% type.typeLists = typeLists;
% type.subTypeLists = subTypeLists;
tis.cells.type = type;

%% cell position

tis.cells.anchors = tis.obI.cell.anchors;

%% insert dat info

for i = 1:length(dat.cid)
   targ = find(tis.cells.cids == dat.cid(i));
   if ~isempty(targ)
      
       if dat.type(i)
           tis.cells.type.typeID(targ) = dat.type(i);
       end
       if dat.subType(i)
           tis.cells.type.subTypeID(targ) = dat.subType(i);
       end
       if sum(dat.onOff(i,:),2)
           tis.cells.type.onOff(targ,:) = dat.onOff(i,:);
       end
       if dat.tracing(i)
           tis.cells.tracing(targ) = dat.tracing(i);
       end
       if dat.seedPosition(i)
           tis.cells.anchors(targ) = dat.seedPosition(i);
       end
       
   end  
    
end


%% Synapses
if ~isempty(obI.nameProps.edges)
tis.syn.edges = obI.nameProps.edges; 
tis.syn.synProp = obI.nameProps.synProp;
tis.syn = getSynMat(tis);
end


%% Format Report

cellHeader = {'CID', 'cell type', 'ON/OFF', 'subtype','tracing','seedPosition'};


rCid = tis.cells.cids';
cellNum = length(rCid);
rTypeID = tis.cells.type.typeID';
rSubTypeID = tis.cells.type.subTypeID';
rOnOff = tis.cells.type.onOff;
rTracing = tis.cells.tracing;
rAnchors = tis.cells.anchors;

clear cTypeNames cCid cSubTypeID cTracing cOnOff cAnchors
for i = 1:cellNum
    
    cCid{i,1} = num2str(rCid(i));
    if rTypeID(i)>0
        cTypeNames{i,1} = tis.cells.type.typeNames{rTypeID(i)};
    else
        cTypeNames{i,1} = 'unassigned';
    end
    
    if (rSubTypeID(i)>0) & (rTypeID(i) > 0)
        cSubTypeID{i,1} = tis.cells.type.subTypeNames{rTypeID(i)}{rSubTypeID(i)};
    else
        cSubTypeID{i,1} = 'unassigned';
    end
        
    if rOnOff(i,1) & rOnOff(i,2)
        cOnOff{i,1} = 'ON\OFF';
    elseif rOnOff(i,1)
        cOnOff{i,1} = 'ON';
    elseif rOnOff(i,2)
        cOnOff{i,1} = 'OFF';
    else
        cOnOff{i,1} = 'unassigned';
    end
    
    cTracing{i,1} = tis.cells.tracing(i);
    cAnchors{i,1} = tis.cells.anchors(i,:);
    

end


cDat = cat(2,cCid,cTypeNames,cSubTypeID,cOnOff,cTracing,cAnchors);
cDat = cat(1,cellHeader,cDat);
   tis.report.cellDat = cDat;


cellSynHeader = {'CID','Cell Type', 'Sub Type', 'Syn In','BipSyn In','Syn Out','RGCSyn out'};


   clear cSyn totSyn
 for i = 1:cellNum
     targ = rCid(i);
    isPre = find(tis.syn.post == targ);
    preClass = tis.syn.preClass(isPre);
    isBip = find(preClass==7);
    
    isPost = find(tis.syn.pre == targ);
    postClass = tis.syn.postClass(isPost);
    isRGC = find(postClass==1);
    
    cSyn(i,:) = {length(isPre),length(isBip),length(isPost),length(isRGC)};
    totSyn(i) = length(isPre) + length(isPost);
 end
 
 [s idx]= sort(totSyn,'descend');
 
 sDat  = cat(2,cCid,cTypeNames,cSubTypeID,cSyn);
 sDat = cat(1,cellSynHeader,sDat(idx,:));
 
 
 tis.report.cellSynDat = sDat;   
 
 
 
 
%% save

save([WPN 'tis.mat'],'tis');




