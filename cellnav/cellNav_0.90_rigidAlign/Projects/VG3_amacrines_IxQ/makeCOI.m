function[] = makeCOI()

global tis glob
datFold = [glob.datDir 'Analysis\Data\preproc\'];
if ~exist(datFold,'dir'),mkdir(datFold);end
clear COI

%% Get data
isAMC = find(tis.cells.type.typeID == 8);
isVGC = intersect(isAMC,find(tis.cells.type.subTypeID == 1));
isTH2 = intersect(isAMC,find(tis.cells.type.subTypeID == 2));
th2Cids = tis.cids(isTH2);
vgcCids = tis.cids(isVGC);
isBip = find(tis.cells.type.typeID == 7);
bipCids = tis.cids(isBip);
bipSub = tis.cells.type.subTypeID(isBip);
bipSubNames = tis.cells.type.subTypeNames{7};
isRGC = find(tis.cells.type.typeID == 1);
rgcCids = tis.cids(isRGC);
rgcSub = tis.cells.type.subTypeID(isRGC);
isUnk = find(tis.cells.type.typeID == 0);
unkCids = [0 tis.cids(isUnk)];

rgcSubNames = tis.cells.type.subTypeNames{1};
bpcSubNames = tis.cells.type.subTypeNames{7};

COI.cids = tis.cids;
COI.isAMC = isAMC;
COI.isVGC = isVGC;
COI.isTH2 = isTH2;
COI.vgcCids = vgcCids;
COI.th2Cids = th2Cids;
COI.isBip = isBip;
COI.bipCids = bipCids;
COI.bipSub = bipSub;
COI.bipSubNames = bipSubNames;
COI.isRGC = isRGC;
COI.rgcSub = rgcSub;
COI.rgcSubNames = rgcSubNames;
COI.bpcSubNames = bpcSubNames;
COI.unkCids = unkCids;





%% Bip subtypes
COI.onBipSubs = [6 7 8 9 10 11 12 14];
COI.offBipSubs = [1 2 3 4 5 15];


clear COI.bpcGroupLabel COI.bpcGroupNames
% 
% for i = 1:length(COI.bpcSubNames)
%     COI.bpcGroupLabel{i} = COI.bpcSubNames{i};
%     COI.bpcGroupNames{i} = COI.bpcSubNames(i);
% end

COI.bpcGroupLabel{1} = 'off';
COI.bpcGroupNames{1} = {'bc1' 'bc2' 'bc3a' 'bc3b' 'bc4' 'boff' 'bcoff' 'off' 'bc3'};
COI.bpcGroupLabel{2} = 'on';
COI.bpcGroupNames{2} = {'bc5i' 'bc5o' 'bc5t' 'xbc' 'bc6' 'bc7' 'bc8' 'rbc' 'bon' 'bcon' 'on' 'bc5'};
COI.bpcGroupLabel{3} = 'transient';
COI.bpcGroupNames{3} = {'bc5i' 'bc5o' 'bc5t' 'xbc' 'bc5' 'bc3a' 'bc3b' 'bc4' 'bc3'};
COI.bpcGroupLabel{4} = 'sustained';
COI.bpcGroupNames{4} ={'bc6' 'bc7' 'bc8' 'rbc' 'bc1' 'bc2' } ;
COI.bpcGroupLabel{5} = 'unk';
COI.bpcGroupNames{5} = {'bcunk' 'unk' ' '};
COI.bpcGroupLabel{6} = 'bc1';
COI.bpcGroupNames{6} = {'bc1'};
COI.bpcGroupLabel{7} = 'bc2';
COI.bpcGroupNames{7} = {'bc2'};
COI.bpcGroupLabel{8} = 'bc3a';
COI.bpcGroupNames{8} = {'bc3a'};
COI.bpcGroupLabel{9} = 'bc3b';
COI.bpcGroupNames{9} = {'bc3b'};
COI.bpcGroupLabel{10} = 'bc3';
COI.bpcGroupNames{10} = {'bc3'};
COI.bpcGroupLabel{11} = 'bc4';
COI.bpcGroupNames{11} = {'bc4'};
COI.bpcGroupLabel{12} = 'bc5i';
COI.bpcGroupNames{12} = {'bc5i'};
COI.bpcGroupLabel{13} = 'bc5o';
COI.bpcGroupNames{13} = {'bc5o'};
COI.bpcGroupLabel{14} = 'bc5t';
COI.bpcGroupNames{14} = {'bc5t'};
COI.bpcGroupLabel{15} = 'xbc';
COI.bpcGroupNames{15} = {'xbc'};
COI.bpcGroupLabel{16} = 'bc5';
COI.bpcGroupNames{16} = {'bc5'};
COI.bpcGroupLabel{17} = 'bc6';
COI.bpcGroupNames{17} = {'bc6'};
COI.bpcGroupLabel{18} = 'bc7';
COI.bpcGroupNames{18} = {'bc7'};
COI.bpcGroupLabel{19} = 'bc8';
COI.bpcGroupNames{19} = {'bc8'};
COI.bpcGroupLabel{20} = 'bc11';
COI.bpcGroupNames{20} = {'bc11'};

COI.checkBPCLabels = {'bc3a' 'bc3b' 'bc4' 'xbc' 'bc5i' 'bc5o' 'bc5t'  'bc6'  };






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

if length(vgcCids) % If there are VGCs
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

end


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
rgcGroupLabel{4} = '25';
rgcGroupNames{4} = {'25'};
rgcGroupLabel{5} = '3i';
rgcGroupNames{5} = {'3i'};
rgcGroupLabel{6} = '4';
rgcGroupNames{6} = {'4'};
rgcGroupLabel{7} = '4i';
rgcGroupNames{7} = {'4i'};
rgcGroupLabel{8} = '4on';
rgcGroupNames{8} = {'4on'};
rgcGroupLabel{9} = '4ow';
rgcGroupNames{9} = { '4ow'};
rgcGroupLabel{10} = '28'; %on off
rgcGroupNames{10} = {'28'};
rgcGroupLabel{11} = '37'; % on off
rgcGroupNames{11} = {'37' '37c' '37d' '37r' '37v'};
rgcGroupLabel{12} = '5ti';
rgcGroupNames{12} = {'5ti'};
rgcGroupLabel{13} = '5si';
rgcGroupNames{13} = {'5si'};
rgcGroupLabel{14} = '5to';
rgcGroupNames{14} = {'5to'};
rgcGroupLabel{15} = '5so';
rgcGroupNames{15} = {'5so'};
rgcGroupLabel{16} = '51'; %of off
rgcGroupNames{16} = {'51'};
rgcGroupLabel{17} = '6t'; %of off
rgcGroupNames{17} = {'6t'};
rgcGroupLabel{18} = '6sw'; %of off
rgcGroupNames{18} = {'6sw'};
rgcGroupLabel{19} = '6sn'; %of off
rgcGroupNames{19} = {'6sn'};
rgcGroupLabel{20} = '63'; %on off
rgcGroupNames{20} = {'63'};
rgcGroupLabel{21} = '7i';
rgcGroupNames{21} = {'7i' '7id' '7ir' '7iv'};
rgcGroupLabel{22} = '7o';
rgcGroupNames{22} = {'7o'};
rgcGroupLabel{23} = '72';%on off
rgcGroupNames{23} = {'72'};
rgcGroupLabel{24} = '85';
rgcGroupNames{24} = {'85'};
rgcGroupLabel{25} = '8w';
rgcGroupNames{25} = {'8w'};
rgcGroupLabel{26} = 'm3';
rgcGroupNames{26} = {'m3'};


COI.rgcGroupNames = rgcGroupNames;
COI.rgcGroupLabel = rgcGroupLabel;
COI.checkRGCLabels =  { '2an' '4i'  '4on'  '4ow' '3i' '37'  '5si' '5to' ...
    '25' '28'  '5so' '5ti' '63' '85'  '7i' '6sw' '6sn' '6t' '8w' };





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








