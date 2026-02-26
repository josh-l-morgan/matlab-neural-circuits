global tis glob

datFold = [glob.datDir 'Analysis\Data\'];


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
    runCids = useVgc;%unique(roiCid);%MOI.cids;
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


%% Bip to RGC influence

%%load SM
smDir = 'E:\IxQ_KarlsRetinaVG3_2019\Analysis\SMs\'
influenceMat = zeros(length(COI.bpcGroupCids), length(COI.rgcGroupCids));
influenceMatEach = zeros(length(COI.bpcGroupCids), length(COI.rgcGroupCids),length(useVgc));
storeInfluenceOnRgcSyn = cell(length(COI.bpcGroupCids),length(COI.rgcGroupCids), length(useVgc));

lc = 16;
for v  = 1:length(useVgc)
        disp(sprintf('checking vgc cid%d. %d of %d.',useVgc(v),v,length(useVgc)))
        sm = sms(v).sm;
        d = sm.syn2Skel.syn2SynDist;
        %d = sm.skel2skel.linDist;

        pre = sm.syn.pre;
        post = sm.syn.post;
        
        %% get all synapse to node distances
       
        W = exp(-d/lc); % Apply length constant
        
        
        %% run groups
        verbose = 0;
        for rG = 1:length(COI.rgcGroupCids)
            rCids = COI.rgcGroupCids{rG};
            
             
            for bG = 1:length(COI.bpcGroupCids)
                bCids = COI.bpcGroupCids{bG};
             trackPost = zeros(1,length(pre),1); 
                
                storeEach = [];
                %%run cells
                for b = 1:length(bCids)
                    bCid = bCids(b);
                    for r = 1:length(rCids)
                        gCid = rCids(r);
                        
                        isPre = find(pre == bCid);
                        isPost = find(post == gCid);
                        
                        if ~isempty(isPre) & ~isempty(isPost)
                            
                            cutEDist = W(isPre,isPost);
                            if verbose
                                rgcGroupLabel{rG}
                                COI.bpcGroupLabel{bG}
                                cutEDist
                                pause
                            end
                            
                            influenceMatEach(bG,rG,v) = influenceMatEach(bG,rG,v) + sum(cutEDist(:));
                           trackPost(isPost) = trackPost(isPost) + sum(cutEDist,1);
       
                            pause(.01)
                        end
                        
                    end %run each rgc cid
                    
                end % run each bip cid
                storeInfluenceOnRgcSyn{bG,rG,v} = trackPost;

            end % run each bpc type
        end %% run each rgc type
         
    
end
influenceMat = sum(influenceMatEach,3);

%% extract influence matrix
clf
useBGroup = find(sum(influenceMat,2));
useBGroup = [3 4 5 6 7 8 9 20];
useRGroup = find(sum(influenceMat,1));
%useRGroup = [ 1 2 3 5  7  8  9 10  13 14 15 17];
COI.bpcGroupLabel(useBGroup)
COI.rgcGroupLabel(useRGroup)

useInfluenceMat = influenceMat(useBGroup,useRGroup');

subplot(4,1,1:2)
maxInfluence = max(useInfluenceMat(:));
colormap(jet(100))
image(useInfluenceMat * 100/maxInfluence);
title('Bipolar (y) to RGC (x) influence normalized for RGC')
normInfluence = useInfluenceMat./repmat(sum(useInfluenceMat,1),...
    [size(useInfluenceMat,1) 1]);
image(normInfluence*100./max(normInfluence(:)))

sumInfluenceRGC = sum(useInfluenceMat,1);
sumInfluenceBPC = sum(useInfluenceMat,2);

subplot(4,1,3)

bar(sumInfluenceRGC,'k','edgeColor','none')
title('total RGC')
subplot(4,1,4)

bar(sumInfluenceBPC,'k','edgecolor','none')
title('total Bipolar')

%% Recheck merges when changes
clf
onInfluence = sum(normInfluence(4:end,:),1);
offInfluence = sum(normInfluence(1:3,:),1);
OFFbias = (offInfluence-onInfluence)./(offInfluence + onInfluence);
bar(OFFbias)

bar(onInfluence);
bar(offInfluence);

onInfluence = sum(useInfluenceMat(4:end,:),1);
offInfluence = sum(useInfluenceMat(1:3,:),1);
OFFbias = (offInfluence)./(offInfluence + onInfluence);
bar(OFFbias)
ylim([0 1])


%% Get influence stats by RGC type

totInfluence = zeros(length(useBGroup),length(useRGroup));
clear sMeanPol sN
clf
rangeP = [0:.1:1];
for r = 1:length(useRGroup)
    onStack = [];
    for b = 1:length(useBGroup)
        allV = [];
        for v = 1:length(useVgc)
            allV = cat(2,allV,storeInfluenceOnRgcSyn{useBGroup(b),useRGroup(r),v});
        end
        onStack(b,:) = allV;
    end
    
    totInfluence(:,r) = sum(onStack,2);
    
    sumB = sum(onStack,1);
    onStack = onStack(:,sumB>0);
    
    normInfluenceS = onStack./repmat(sum(onStack,1), [size(onStack,1) 1]);
    
    onInfluenceS = sum(normInfluenceS(4:end,:),1);
    offInfluenceS = sum(normInfluenceS(1:3,:),1);
    OFFbiasS = (offInfluenceS)./(offInfluenceS + onInfluenceS);
    sMeanPol(r) = mean(OFFbiasS);
    sN(r) = length(OFFbiasS);
   
    
    histP = hist(OFFbiasS,rangeP);
    subplot(length(useRGroup),1,r)
    bar(rangeP,histP)
    axis on
    ylim([0 40])
    
end

%% Stratification









