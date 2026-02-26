


load('MPN.mat')
load('WPN.mat')

load([MPN 'obI.mat'])
load([MPN 'dsObj.mat'])
load([WPN 'tis.mat'])



%% Get VG3

amcList = find(tis.cells.type.typeID == 8); % list Amacrines
vgcList = tis.cells.type.subTypeLists.vgc;


%% inf: cid, tisNum, cell type, subtype, length, input (total, by type), output (total, by type)
tableName = {'tisNum','cid','cell type', 'sub type', 'length', 'all in', 'bip in',...
    'amc in', 'unk in', 'all out', 'bip out', 'amc out', 'rgc out', 'unk out'}

cellNum = length(tis.cells.cids);
tis.cells.type.typeNames
tableDat = zeros(cellNum,length(tableName));
for i = 1:cellNum
    
    cid = tis.cells.cids(i);
    tableDat(i,1) = i;
    tableDat(i,2) = cid;
    tableDat(i,3) = tis.cells.type.typeID(i);
    tableDat(i,4) = tis.cells.type.subTypeID(i);
    tableDat(i,5) = 0;
    
    isPre = tis.syn.pre == cid;
    isPost = tis.syn.post == cid;
    
    tableDat(i,6) = sum(isPost);
    preTypes = tis.syn.preClass(isPost);
    tableDat(i,7) = sum(preTypes == 7);
    tableDat(i,8) = sum(preTypes == 8);
    tableDat(i,9) = sum(preTypes == 0);
    
    tableDat(i,10) = sum(isPre);
    postTypes = tis.syn.postClass(isPre);
    tableDat(i,11) = sum(postTypes == 7);
    tableDat(i,12) = sum(postTypes == 8);
    tableDat(i,13) = sum(postTypes == 1);
    tableDat(i,14) = sum(postTypes == 0);
    
end

%% Sort table
% useDat = (tableDat(:,3) == 8) & ( tableDat(:,4) == 1);
% [sortDat ord] = sort(tableDat(:,2));
% useOrd = ord(useDat);

cidOrd = [2 3 4 10 11 13 14 5];
useOrd = [];
for i = 1:length(cidOrd)
   useOrd = [useOrd tableDat(tableDat(:,2)==cidOrd(i),1)];
end

useTab = tableDat(useOrd,:);

%% Get cell stats
for i = 1:size(useTab,1)
   
    cid = useTab(i,2);
    smFile = sprintf('%s\\SMs\\sm_cid%d.mat',WPN,cid);
    if 0%~exist(smFile,'file')
        makeSM(cid);
    end
    
    
        load(smFile)
        
        if 0
            clf
            renderFV(sm.nep.fv,[1 0 0],.2)
            %showRadSurf(nep.pos,nep.edges, nep.nodeRad,[0 1 0],.4)
            showRadSurf(sm.nep.pos,sm.nep.edges, sm.nep.meanNodeRad, [ 0 0 1],.2)
            showRadSurf(sm.nep.pos,sm.nep.edges, sm.nep.meanNodeRad * 0 + .05, [ 0 1 0],1)
                        
        end
    
    useTab(i,5) = sum(sm.arbor.edgeProps.length);
    
end


%% convert to cell
cellTab = num2cell(useTab);%size(useTab,1),size(useTab,2))
for i = 1:size(cellTab,1)
    if useTab(i,3) == 0
        cellTab{i,3} = '?';
    else
        cellTab{i,3} = tis.cells.type.typeNames{useTab(i,3)};
    end
    if useTab(i,4) ==0;
        cellTab{i,4} = '?';
    else
        cellTab{i,4} = tis.cells.type.subTypeNames{useTab(i,4)};
    end
end

fullTab = cat(1,tableName,cellTab)

%{

tis.cells.type.subTypeNames{
%}


%% Get cell stats














