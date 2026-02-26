%%Collect synapse number and arbor size information from each vGlut3 cell
%%for comparison

compType = 2; %3 , 2
%%3 = RGCs


clf

%% load data
if 0

    global glob tis

    figure
    SPN = [glob.datDir 'Analysis\Data\preproc\'];
    load([SPN 'ptDat.mat']);
    load([SPN 'ROI.mat']);
    %load([SPN 'ROI2.mat']);
    load([SPN 'SOI.mat']);
    load([SPN 'GOI.mat']);
    load([SPN 'NOI.mat']);
    load([SPN 'MOI.mat']);
    load([SPN 'COI.mat']);



    %%Load all sms
    smDir = [glob.dir.Volumes  glob.vol.activeName '\Analysis\SMs\'];
    clear sms
    roiCid = ptDat(:,3);
    runCids = unique(roiCid);%MOI.cids;
    for i = 1:length(runCids);
        cid = runCids(i);
        %%Get distances between nodes
        disp(sprintf('loading data for cell %d.  Cell %d of %d.',cid,i,length(runCids)));
        fileName = sprintf('sm_cid%d.mat',cid);
        useSM(i) = 1;
        if exist([smDir fileName],'file')
            sm = load([smDir fileName]);
            load([smDir fileName]);
            sm.nep.swcS = nep2swc(sm.nep);
            sms(i).sm = sm;
        end
    end
    useVgc = runCids(useSM>0);

    s = makeSynapseClassifyer(COI); %make structure describing types of synapses
end

useVgc =[ 2 3 4 5 13 14];

%% set variables
distBin = 3;
bRange = [0 : .1 : 20];
ignoreSelf = 0; %ignore synapses between the same sets of cell IDs, 2 = only self


%% colect data on cells

cellSynCount = zeros(length(s),length(useVgc));
cellCount = cellSynCount;
clear allPart;
countG = zeros(length(s),1);
for v = 1:length(useVgc)

    cid = useVgc(v);
    vTarg = find(runCids == cid);
    useSM(v) = 1;
    disp(sprintf('loading sm %d, %d of %d',cid,v,length(useVgc)))
    sm = sms(vTarg).sm;%load([smDir fileName],'syn','syn2Skel','nep');

    syn = sm.syn; % get syn information
    synD = sm.syn2Skel.syn2SynDist;
    skelD = sm.syn2Skel.syn2SkelDist;
    nodeLengths = sm.nep.props.nodeLength;

    %% Create synapse groups filled with relevant IDs using COI and subtypes


    %%Collect synapses
    for i = 1:length(s)
        hit = [1:length(syn.pre)]';
        
        if ~isempty(s(i).synType)
            hit = intersect(hit,find(syn.synType == s(i).synType));
        else
            hit = intersect(hit,find(syn.synType<inf)) ;
        end

        if isempty(s(i).input)
            part = syn.pre(hit);

        elseif s(i).input
            hit = intersect(hit, find(syn.post == cid));
            part = syn.pre(hit);

        else
            hit = intersect(hit, find(syn.pre == cid));
            part = syn.post(hit);
        end


        if ~s(i).checkCids
            s(i).synIds = hit;

        else
            isPart = part * 0;
            for h = 1:length(hit)
                if sum(s(i).cids == part(h))
                    isPart(h) = 1;
                end
            end
            s(i).synIds = hit(isPart>0);
            part = part(isPart'>0);
        end
        cellCount(i,v) = length(unique(part));
        allPart{i,v} = part;
        cellSynCount(i,v) = length(s(i).synIds);
    end


end

%% Combine IDs from different vg3s
clear countCids
for i = 1:size(allPart,1)
    iPart = [];
    for v = 1:size(allPart,2)
        iPart = cat(1,iPart,allPart{i,v});
    end
    countCids{i} = iPart;
end


%% choose synTypes to compare

sNames = {s(:).name}'

sNames = {s(:).name}'
clear sg
if compType == 1
    sg(1).names = {'bc3a'};
    sg(2).names = {'bc3b'};
    sg(3).names = {'bc4'};
    checkCol = [0 0 1; 0 1 1; 0 1 1; 1 0 0; 1 .5 0; 1 1 0; 0 0 0];

elseif compType ==2

    checkBPCLabels = COI.checkBPCLabels; 
    for i = 1:length(checkBPCLabels)
        sg(i).names = checkBPCLabels(i);
    end
    
    checkCol = jet(length(sg)) *.6 + .4;
    checkCol = jet(length(sg)+3) *.9 + .1;
    checkCol = checkCol([2 4 5 7 8 9 10 11],:);
    checkCol(checkCol>1) = 1;
    % 
    % 
    % checkCol = hsv(100) *.8 + .2;
    % checkCol = checkCol([70 60 52 21 14 8 1 95],:);
    % checkCol(checkCol>1) = 1;


elseif compType ==3
    
    checkRGCLabels = COI.checkRGCLabels; 
    for i = 1:length(checkRGCLabels)
        sg(i).names = checkRGCLabels(i);
    end
    
    checkCol = jet(length(sg)) *.6 + .4;
    checkCol(checkCol>1) = 1;
% 
%     checkCol = jet(length(sg));
%     dimGroups = [];
%     checkCol(dimGroups,:) = 1- checkCol(dimGroups,:) * .1;

      %checkCol = [0 0 1; 0 .4 1; 0 .8 1; 1 .8 0; 1 .6 0; ,1 .4 0; ,1 .2 0 ; .8 0 .1; 1 0 0];

elseif compType == 4 %% basic syn types

 
    sg(2).names = {'non-rib'};
    sg(1).names = {'rib'};
    sg(3).names = {'toAMCs' 'toUnk'};
    sg(4).names = {'RGCs'};
    checkCol = [.6 1 .6; 1 .6 .6; 1 .6 1; .6 .6 1];

elseif compType == 5

    sg(1).names = {'bc3a' 'bc3b' 'bc4'};
    sg(2).names = {'bc5i' 'bc5o' 'bc5t' 'xbc' 'bc6' 'bc5i'};
 
elseif compType == 6
   sg(1).names = {'1wt' '2aw' '2an' '4' '4i' '28' '37 37c 37d 37r 37v' '5 5si 5so' '51' ...
       '7 7id 7ir 7iv 7o' '8 8w 8n'};
    sg(2).names = {'4ow'};
    sg(3).names = {'63'};
   

    checkCol = [.2 .2 .2; .7 .3 .1; .1 .3 .9];
elseif compType == 7
    sg(1).names = {'RGCs'};
    sg(2).names = {'toAMCs'};
    sg(3).names = {'allOutputs'};


    checkCol = [.2 .2 .2; .7 .3 .1; .1 .3 .9];

elseif compType == 8 %All RGC types

    useRGCs = [21:length(sNames)]
    for u = 1:length(useRGCs)
            sg(u).names = sNames(useRGCs(u));
    end
    checkCol = hsv(length(useRGCs));

elseif compType == 9 %Pool identified RGCs
     sg(1).names = {'1wt' '2aw' '2sn' '4' '4i' '4ow' '28' '37 37c 37d 37r 37v'...
     '5 5si 5so' '51' '63' '7 7id 7ir 7iv 7o' '8 8w 8n'};
     checkCol = [1 0 0];

elseif compType == 10
      
    sg(1).names = {'bc1'};
    sg(2).names = {'bc2'};
    sg(3).names = {'bc3a'};
    sg(4).names = {'bc3b'};
    sg(5).names = {'bc4'};
    sg(6).names = {'xbc'};
    sg(7).names = {'bc5i'};
    sg(8).names = {'bc5o'};
    sg(9).names = {'bc5t'};
    sg(10).names = {'bc6'};
    sg(11).names = {'bc7'};
    sg(12).names = {'bc8/9'};
    
    checkCol = hsv(length(sg));
       checkCol = jet(length(sg));
    dimGroups = [];
    checkCol(dimGroups,:) = 1- checkCol(dimGroups,:) * .1;
   

end

%% pick s's for sgs
expNames = cat(1,{sg(:).names});
for g = 1:length(sg)
    sg(g).sT = [];
    nams = sg(g).names;
    for n = 1:length(nams)
        m1 = find(strcmp(sNames,sg(g).names{n}));
        if ~isempty(m1)
            sg(g).sT(n) = m1;
        else
            sg(g).sT(n) = 0;
        end
    end
end

%%Collect data from cellSynCount
cellGroupCount = zeros(length(sg),length(useVgc));
cellGroupCidCount = cellGroupCount;
clear gCidNum
for g = 1:length(sg)
    if sum(sg(g).sT)
        cellGroupCount(g,:) = sum(cellSynCount(sg(g).sT,:),1);
        cellGroupCidCount(g,:) = sum(cellCount(sg(g).sT,:),1);
        gCids = [];
        for t = 1:length(sg(g).sT)
            gCids = cat(1,gCids,countCids{sg(g).sT(t)});
        end
        gCidNum(g) = length(unique(gCids));
    else
        gCidNum(g) = 0;
    end
end

allGroupCount = sum(cellGroupCount,2);
isHit = find(allGroupCount>0);
allGroupCidCount = sum(cellGroupCidCount,1);

%% Show pie charts
clf
for i = 1:length(useVgc)
    subplot(2,3,i)
    vals = cellGroupCount(isHit,i);
    titleStr = sprintf('cid = %d, n = %d',useVgc(i),sum(vals))
    pChart = pie(vals)
    for p = 1:length(vals)
        %pChart(p*2-1).FaceColor
        pChart(p*2-1).FaceColor = checkCol(p,:);
    end
        title(titleStr)
    if i == length(useVgc)
        clear leg
        for n = 1:length(isHit)
            nam = expNames{isHit(n)};
            toe = [];
            for c = 1:length(nam)
                toe = [toe ' ' nam{c}];
            end

            leg{n} = toe;
        end
        legend(leg)
    end
end




%% Show combined

clf
    vals = allGroupCount(isHit);
    titleStr = sprintf('allCells, n = %d',sum(vals))
    pChart = pie(vals);
    for p = 1:length(vals)
        %pChart(p*2-1).FaceColor
        pChart(p*2-1).FaceColor = checkCol(p,:);
    end
        title(titleStr)
        clear leg
        for n = 1:length(isHit)
            nam = expNames{isHit(n)};
            toe = [];
            for c = 1:length(nam)
                toe = [toe ' ' nam{c}];
            end

            leg{n} = toe;
        end
        legend(leg)

        for i = 1:length(allGroupCount)
            fprintf('%s    %d \n',expNames{i}{1},allGroupCount(i))
        end

        return

%%
clf
useVgc
%useCell = [1 2 3 4  7 8];
useCell = 1:length(useVgc);

sampCount = cellGroupCount(:,useCell);
useCount = [1:size(sampCount,1)];


allSampCount = sum(sampCount,2);
a = sampCount./sum(sampCount,1)*100;
ma = mean(a,2);
sea = std(a')/sqrt(length(useCell));



gNames = {sg(:).names }



sumCount = sum(sampCount,1);
sampCount = sampCount./repmat(sumCount,[size(sampCount,1) 1])*100;
res = bootMedInt(sampCount)
ma = res(:,2);

meanCount = mean(sampCount,2);
seCount = std(sampCount,[],2)/sqrt(size(sampCount,2));

hold on
seW = .1;
mW = .3;
for i = 1:size(cellGroupCount,1);
    
    line([i i],[meanCount(i)-seCount(i) meanCount(i)+seCount(i)],...
        'color','k','linewidth',2)
    line([i-seW i+seW],[meanCount(i) + seCount(i) meanCount(i) + seCount(i)],...
        'color','k','linewidth',2)
    line([i-seW i+seW],[meanCount(i) - seCount(i) meanCount(i) - seCount(i)],...
        'color','k','linewidth',2)
    line([i-mW i+mW],[meanCount(i)  meanCount(i) ],...
        'color','k','linewidth',2)
    
    sc = swarmchart(ones(length(useCell),1) * i,sampCount(i,:),...
        'sizedata',100,'markerfacecolor',checkCol(i,:),'markerfacealpha',.5,...
        'xjitterwidth',.3,'xjitter','density','linewidth',2,'markeredgecolor','k');
    %legend(leg)
    xlim([0 max(useCount)+1])
end



%% Print figure
if 0
    %fDir = uigetdir;
    filename = [fDir '\' 'RGCfractions_2022update']
    set(gcf,'renderer','Painters')
    print('-depsc','-tiff','-r300', '-painters',[filename,'.eps'])
    
end







