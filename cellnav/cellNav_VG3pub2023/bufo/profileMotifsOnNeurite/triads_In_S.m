%%Characterize triads using s structure for synapses


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
        %sm = load([smDir fileName],'skel2skel','nep','syn2skel','syn');
        load([smDir fileName]);
        sm.nep.swcS = nep2swc(sm.nep);
        sms(i).sm = sm;
    end
    useVgc = runCids(useSM>0);

    s = makeSynapseClassifyer(COI); %make structure describing types of synapses

end

%% cell to cell triads
syn = tis.syn;

vgcCids = COI.vgcCids;

c2cTriads = [];
for i = 1:length(vgcCids)
    synToVGC = find(syn.post == vgcCids(i));
    synFromVGC = find(syn.pre == vgcCids(i));
    cidsToVGC = setdiff(unique(syn.pre(synToVGC)),0);
    cidsFromVGC = setdiff(unique(syn.post(synFromVGC)),0);

    for y = 1:length(cidsToVGC)
        for x = 1:length(cidsFromVGC)
            
            hits =  sum((syn.pre == cidsToVGC(y)) & (syn.post == cidsFromVGC(x)));
            if hits
            c2cTriads = cat(1,c2cTriads,[cidsToVGC(y) vgcCids(i) cidsFromVGC(x)]);
            end
        end
    end
end


autapse = [(c2cTriads(:,1)==c2cTriads(:,2)) | (c2cTriads(:,2)==c2cTriads(:,3)) | (c2cTriads(:,1)==c2cTriads(:,3)) ]; 
c2cTriads = c2cTriads(autapse==0,:);


%%Get triad numbers, order 1 to 2, 1 to 3, 2 to 3
c2cCount = c2cTriads * 0;
for t = 1:size(c2cTriads,1)
    
    s1to2 = sum(syn.pre==c2cTriads(t,1) & (syn.post==c2cTriads(t,2)));
    s1to3 = sum(syn.pre==c2cTriads(t,1) & (syn.post==c2cTriads(t,3)));
    s2to3 = sum(syn.pre==c2cTriads(t,2) & (syn.post==c2cTriads(t,3)));
    c2cCount(t,:) = [s1to2 s1to3 s2to3];
end

%%Get syn to syn distances

%%Get triad types



%% Type to type triads

%checkRGCLabels = {'1wt' '2aw' '2an' '4i' '4ow' '28' '37' '5si' '5ti' '6sw' '6sn' '63' '85' '8w'};
checkRGCLabels = COI.checkRGCLabels; 
checkBPCLabels = COI.checkBPCLabels; 



clear useR useB
%%for each COIgroup determine if it is a label to check
for i = 1:length(checkRGCLabels)
   useR(i) = find(strcmp(COI.rgcGroupLabel,checkRGCLabels{i}));
end
for i = 1:length(checkBPCLabels)
    useB(i) = find(strcmp(COI.bpcGroupLabel,checkBPCLabels{i}));
end

useRNames = {COI.rgcGroupLabel{useR}};
useBNames = {COI.bpcGroupLabel{useB}};
rbMat = zeros(length(useB),length(useR));
clear memRCids memBCids
for r = 1:length(useR); 
    

    rCids = COI.rgcGroupCids{useR(r)};
    memRCids{r} = rCids;
    for b = 1:length(useB)
        bCids = COI.bpcGroupCids{useB(b)};
        if r == 1
            memBCids{b} = bCids;
        end

        for y = 1:length(rCids)
            for x = 1:length(bCids)

                hits = sum((syn.pre==bCids(x)) & (syn.post==rCids(y)));
                rbMat(b,r) = rbMat(b,r) + hits;
            end
        end

        
    end
end

cmap = jet(1000);
cmap(1,:) = 0;
colormap(cmap)
image(rbMat*1000/max(rbMat(:)))
maxTriadNumber = max(rbMat(:))

%% Check arbors
global glob
MPN = glob.NA.MPN;
load([MPN 'dsObj.mat'])
load([MPN 'obI.mat'])


zThresh = 125;


voxCellNum = zeros(length(useR),1);
for r = 1:length(useR)
    useRNames{r}
    clf 
    hold on
    for c = 1:length(memRCids{r})
        cellTarg = find(obI.cell.name == memRCids{r}(c));
        obIds = obI.cell.obIDs{cellTarg};
        subs = cat(1,dsObj(obIds).subs);
        %scatter(subs(:,1),subs(:,3))
        voxCellNum(r,c) = sum(subs(:,3)>zThresh);
    end

end
hold off


voxNum = sum(voxCellNum,2);
voxVol = obI.em.dsRes^3;
typeVol = voxNum.*voxVol;
rbMatS = rbMat./ repmat(typeVol',[size(rbMat,1) 1]);
rbMatS(isnan(rbMatS)) = 0;

zThreshBPC = 350;
voxCellNumBPC = zeros(length(useB),1);
for r = 1:length(useB)
    for c = 1:length(memBCids{r})
        cellTarg = find(obI.cell.name == memBCids{r}(c));
        obIds = obI.cell.obIDs{cellTarg};
        subs = cat(1,dsObj(obIds).subs);
        %scatter(subs(:,1),subs(:,3)), 
        if ~isempty(subs)
        voxCellNumBPC(r,c) = sum(subs(:,3)<zThreshBPC);
        else
            disp(sprintf('cid %d is empty but has synapse',memBCids{r}(c)))
        end
    end
end
hold off



voxNumBPC = sum(voxCellNumBPC,2);
typeVolBPC = voxNumBPC.*voxVol;

subplot(4,1,1:2)


image(rbMatS*1000/max(rbMatS(:)))
image(rbMat*1000/max(rbMat(:)))

xticks([1:length(checkRGCLabels)])
xticklabels(checkRGCLabels)

yticks([1:length(checkBPCLabels)])
yticklabels(checkBPCLabels)



subplot(4,1,3)
bar(typeVol,'k','edgeColor','none')



subplot(4,1,4)
bar(typeVolBPC,'k','edgeColor','none')










