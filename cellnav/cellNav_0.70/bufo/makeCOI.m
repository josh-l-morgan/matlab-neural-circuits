function[] = makeCOI()

global tis glob

datFold = [glob.datDir 'Analysis\Data\'];
if ~exist(datFold,'dir'),mkdir(datFold);end

%% Get data
isAMC = find(tis.cells.type.typeID == 8); 
isVGC = intersect(isAMC,find(tis.cells.type.subTypeID == 1));
vgcCids = tis.cids(isVGC);
isBip = find(tis.cells.type.typeID == 7); 
bipCids = tis.cids(isBip);
bipSub = tis.cells.type.subTypeID(isBip);
bipSubNames = tis.cells.type.subTypeNames{7};
isRGC = find(tis.cells.type.typeID == 1);
rgcCids = tis.cids(isRGC);
rgcSub = tis.cells.type.subTypeID(isRGC);

rgcSubNames = tis.cells.type.subTypeNames{1};
bpcSubNames = tis.cells.type.subTypeNames{7};


COI.isAMC = isAMC;
COI.isVGC = isAMC;
COI.vgcCids = vgcCids;
COI.isBip = isBip;
COI.bipCids = bipCids;
COI.bipSub = bipSub;
COI.bipSubNames = bipSubNames;
COI.isRGC = isRGC;
COI.rgcSub = rgcSub;
COI.rgcSubNames = rgcSubNames;
COI.bpcSubNames = bpcSubNames;


%% Bip subtypes
COI.onBipSubs = [6 7 8 9 10 11 12 14];
COI.offBipSubs = [1 2 3 4 5 15];


clear COI.bpcGroupLabel COI.bpcGroupNames
for i = 1:length(COI.bpcSubNames)
    COI.bpcGroupLabel{i} = COI.bpcSubNames{i};
    COI.bpcGroupNames{i} = COI.bpcSubNames(i);
end

clear COI.bpcGroupID
for i = 1: length(COI.bpcGroupLabel);
    COI.bpcGroupID{i} = [];
    for r = 1:length(COI.bpcGroupNames{i})
        for n = 1:length(COI.bpcSubNames)
            if strcmp(COI.bpcGroupNames{i}{r},COI.bpcSubNames{n})
                COI.bpcGroupID{i} =  cat(2,COI.bpcGroupID{i},n);
                break
            end
        end
    end
end


COI.synToVGC = [];
for i = 1:length(vgcCids)
    COI.synToVGC = cat(1,COI.synToVGC,find(tis.syn.post == COI.vgcCids(i)));
end

COI.preToVGC = tis.syn.pre(COI.synToVGC);
COI.typesPreToVGC = tis.syn.preClass(COI.synToVGC);
COI.cidsPreToVGC = unique(COI.preToVGC);

clear COI.cPreType COI.sPreType
for i = 1:length(COI.cidsPreToVGC);
    targ = find(tis.cids == COI.cidsPreToVGC(i));
    if isempty(targ)
        COI.cPreType(i) = 0;
        COI.sPreType(i) = 0;
    else
       COI.cPreType(i) = tis.cells.type.typeID(targ); 
       COI.sPreType(i) = tis.cells.type.subTypeID(targ); 
    end
end

COI.bpcPreToVGC = COI.cidsPreToVGC(COI.cPreType == 7);
COI.bpcSubs = COI.sPreType(COI.cPreType == 7);


clear COI.bpcGroupCids
for g = 1:length(COI.bpcGroupID)
    isType = [];
    bCids = [];
    for i = 1:length(COI.bpcGroupID{g})
       hit =  find((COI.bpcSubs == COI.bpcGroupID{g}(i)));
       bCids = cat(1,bCids,COI.bpcPreToVGC(hit));
    end
     COI.bpcGroupCids{g} = bCids;
end
        
isKnown = find((COI.bpcSubs>0)&(COI.bpcSubs<=12));
isNotKnown = find((COI.bpcSubs==0)|(COI.bpcSubs>12));
COI.identifiedBips = unique(COI.bpcPreToVGC(isKnown));
COI.unidentifiedBips = unique(COI.bpcPreToVGC(isNotKnown));

COI.knownBipSubs = tis.cells.type.subTypeNames{7}(COI.bpcSubs(isKnown));


useSubs = [];
totSubCids = 0;
COI.groupedCids = [];
clear COI.bpcGroupCellSynNum COI.bpcGroupSynNum COI.bpcGroupCellNum
for i = 1:length(COI.bpcGroupCids)
    checkCids = COI.bpcGroupCids{i};
    totSubCids = totSubCids + length(checkCids);
    COI.groupedCids = cat(1,COI.groupedCids,COI.bpcGroupCids{i});
    
   if ~isempty(checkCids);
       useSubs = [useSubs i];
   end 
   
   COI.bpcGroupCellNum(i) = length(COI.bpcGroupCids{i});
   COI.bpcGroupSynNum(i) = 0;
   for c = 1:length(checkCids)
      COI.bpcGroupCellSynNum{i}(c) =  sum(COI.preToVGC == checkCids(c));
      COI.bpcGroupSynNum(i) = sum(COI.bpcGroupCellSynNum{i});
   end
   
end
COI.ungroupedBPCCids = setdiff(COI.bpcPreToVGC,COI.groupedCids);

%%Make bar
barBpcGroupCells = [];
barBpcGroupSyns = [];
for i = 1:length(COI.bpcGroupCids)
    barBpcGroupCells(length(barBpcGroupCells)+1: ...
        length(barBpcGroupCells) + length(COI.bpcGroupCids{i})) = i;
    barBpcGroupSyns(length(barBpcGroupSyns) + 1: ...
        length(barBpcGroupSyns) + COI.bpcGroupSynNum(i)) = i;
end

colormap(colorcube(length(COI.bpcGroupCids)*2))
% image(barBpcGroupCells')
% image(barBpcGroupSyns')

report = cat(2,COI.bpcGroupLabel',mat2cell(COI.bpcGroupCellNum',ones(length(COI.bpcGroupLabel),1)),...
    mat2cell(COI.bpcGroupSynNum',ones(length(COI.bpcGroupLabel),1)));








%% RGC subtypes
clear rgcGroupLabel
rgcGroupLabel{1} = '1wt';
rgcGroupNames{1} = {'1wt'};
rgcGroupLabel{2} = '2aw';
rgcGroupNames{2} = {'2aw'};
rgcGroupLabel{3} = '2an';
rgcGroupNames{3} = {'2an'};
rgcGroupLabel{4} = '4';
rgcGroupNames{4} = {'4' '4i' '4on' '4ow'};
rgcGroupLabel{5} = '4i';
rgcGroupNames{5} = {'4i'};
rgcGroupLabel{6} = '4on';
rgcGroupNames{6} = {'4on'};
rgcGroupLabel{7} = '4ow';
rgcGroupNames{7} = { '4ow'};
rgcGroupLabel{8} = '28'; %on off
rgcGroupNames{8} = {'28'};
rgcGroupLabel{9} = '37'; % on off
rgcGroupNames{9} = {'37' '37c' '37d' '37r' '37v'};
rgcGroupLabel{10} = '5s';
rgcGroupNames{10} = {'5' '5si' '5so'};
rgcGroupLabel{11} = '5ti';
rgcGroupNames{11} = {'5ti'};
rgcGroupLabel{12} = '5to';
rgcGroupNames{12} = {'5to'};
rgcGroupLabel{13} = '51'; %of off
rgcGroupNames{13} = {'51'};
rgcGroupLabel{14} = '63'; %on off
rgcGroupNames{14} = {'63'};rgcGroupLabel{15} = '7';
rgcGroupNames{15} = {'7' '7id' '7ir' '7iv' '7o'};
rgcGroupLabel{16} = '72';%on off
rgcGroupNames{16} = {'72'}; 
rgcGroupLabel{17} = '8';
rgcGroupNames{17} = {'8' '8w' '8n'};


COI.rgcGroupNames = rgcGroupNames;
COI.rgcGroupLabel = rgcGroupLabel;





clear COI.rgcGroupID
for i = 1: length(rgcGroupLabel);
    rgcGroupID{i} = [];
    for r = 1:length(rgcGroupNames{i})
        for n = 1:length(rgcSubNames)
            if strcmp(rgcGroupNames{i}{r},rgcSubNames{n})
                rgcGroupID{i} =  cat(2,rgcGroupID{i},n);
                break
            end
        end
    end
end
COI.rgcGroupID = rgcGroupID;

COI.synFromVGC = [];
for i = 1:length(vgcCids)
    COI.synFromVGC = cat(1,COI.synFromVGC,find(tis.syn.pre == COI.vgcCids(i)));
end

COI.postToVGC = tis.syn.post(COI.synFromVGC);
COI.typesPostToVGC = tis.syn.postClass(COI.synFromVGC);
COI.cidsPostToVGC = unique(COI.postToVGC);

clear COI.cType COI.sType
for i = 1:length(COI.cidsPostToVGC);
    targ = find(tis.cids == COI.cidsPostToVGC(i));
    if isempty(targ)
        COI.cType(i) = 0;
        COI.sType(i) = 0;
    else
       COI.cType(i) = tis.cells.type.typeID(targ); 
       COI.sType(i) = tis.cells.type.subTypeID(targ); 
    end
end

COI.rgcsPostToVGC = COI.cidsPostToVGC(COI.cType == 1);
COI.rgcSubs = COI.sType(COI.cType == 1);


clear COI.rgcGroupCids
for g = 1:length(COI.rgcGroupID)
    isType = [];
    gCids = [];
    for i = 1:length(rgcGroupID{g})
       hit =  find((COI.rgcSubs == COI.rgcGroupID{g}(i)));
       gCids = cat(1,gCids,COI.rgcsPostToVGC(hit));
    end
     COI.rgcGroupCids{g} = gCids;
end
        
isKnown = find(COI.rgcSubs>4);
knownSubs = tis.cells.type.subTypeNames{1}(COI.rgcSubs(isKnown));


useSubs = [];
totSubCids = 0;
COI.groupedCids = [];
clear groupCellSynNum groupSynNum groupCellNum
for i = 1:length(COI.rgcGroupCids)
    checkCids = COI.rgcGroupCids{i};
    totSubCids = totSubCids + length(checkCids);
    COI.groupedCids = cat(1,COI.groupedCids,COI.rgcGroupCids{i});
    
   if ~isempty(checkCids);
       useSubs = [useSubs i];
   end 
   
   COI.groupCellNum(i) = length(COI.rgcGroupCids{i});
   COI.groupSynNum(i) = 0;
   for c = 1:length(checkCids)
      COI.groupCellSynNum{i}(c) =  sum(COI.postToVGC == checkCids(c));
      COI.groupSynNum(i) = sum(COI.groupCellSynNum{i});
   end
   
end
COI.ungroupedCids = setdiff(COI.rgcsPostToVGC,COI.groupedCids);

%%Make bar
barGroupCells = [];
barGroupSyns = [];
for i = 1:length(COI.rgcGroupCids)
    barGroupCells(length(barGroupCells)+1: ...
        length(barGroupCells) + length(COI.rgcGroupCids{i})) = i;
    barGroupSyns(length(barGroupSyns) + 1: ...
        length(barGroupSyns) + COI.groupSynNum(i)) = i;
end

colormap(colorcube(length(COI.rgcGroupCids)*2))
% image(barGroupCells')
% image(barGroupSyns')

COI.report = cat(2,COI.rgcGroupLabel',mat2cell(COI.groupCellNum',ones(length(COI.rgcGroupLabel),1)),...
    mat2cell(COI.groupSynNum',ones(length(COI.rgcGroupLabel),1)));

save([datFold 'COI.mat'],'COI')








