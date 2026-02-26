global tis glob

datFold = [glob.datDir 'Analysis\Data\'];
load([datFold 'COI.mat']);




%% Bip to RGC influence

%%load SM
smDir = 'E:\IxQ_KarlsRetinaVG3_2019\Analysis\SMs\'
influenceMat = zeros(length(COI.bpcGroupCids), length(COI.rgcGroupCids));
influenceMatEach = zeros(length(COI.bpcGroupCids), length(COI.rgcGroupCids),length(vgcCids));
storeInfluenceOnRgcSyn = cell(length(COI.bpcGroupCids),length(COI.rgcGroupCids), length(vgcCids));

lc = 20;
for v  = 1:length(vgcCids)
    fileName = sprintf('sm_cid%d.mat',vgcCids(v));
    if exist([smDir fileName],'file')
        useSM(v) = 1;
        load([smDir fileName])
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
                sum(trackPost)
                storeInfluenceOnRgcSyn{bG,rG,v} = trackPost;

            end % run each bpc type
        end %% run each rgc type
         
    else
        useSM(v) = 0;
    end
    
    
    
end
influenceMat = sum(influenceMatEach,3);

%% extract influence matrix
useBGroup = find(sum(influenceMat,2));
useBGroup = [3 4 5 6 7 8 9];
useRGroup = find(sum(influenceMat,1));
COI.bpcGroupLabel(useBGroup)
rgcGroupLabel(useRGroup)

useInfluenceMat = influenceMat(useBGroup,useRGroup');

maxInfluence = max(useInfluenceMat(:));
colormap(jet(100))
image(useInfluenceMat * 100/maxInfluence);

normInfluence = useInfluenceMat./repmat(sum(useInfluenceMat,1),...
    [size(useInfluenceMat,1) 1]);
image(normInfluence*100./max(normInfluence(:)))

sumInfluenceRGC = sum(useInfluenceMat,1);
sumInfluenceBPC = sum(useInfluenceMat,2);


bar(sumInfluenceRGC)
bar(sumInfluenceBPC)

%% Recheck merges when changes

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


%% Get influence stats by RGC type

totInfluence = zeros(length(useBGroup),length(useRGroup));
clear sMeanPol sN
clf
rangeP = [0:.1:1];
for r = 1:length(useRGroup)
    onStack = [];
    for b = 1:length(useBGroup)
        allV = [];
        for v = 1:length(vgcCids)
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









