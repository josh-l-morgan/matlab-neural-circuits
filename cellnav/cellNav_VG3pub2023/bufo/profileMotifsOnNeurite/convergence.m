%% Plot number of synapses between pairs of cells.

global glob tis

clf

%% load data
if 0
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


%% identify synapses in groups

syn = tis.syn;
vgcCids = COI.vgcCids;

% %%Collect synapses
% for i = 1:length(s)
%     if ~isempty(s(i).synType)
%         hit = find(syn.synType == s(i).synType);
%     else
%         hit = find(syn.synType<inf) ;
%     end
% 
%     if s(i).input
%         isVsyn = [];
%         for v = 1:length(vgcCids)
%             isVsyn = [isVsyn; find(syn.post == vgcCids(v))];
%         end
%         hit = intersect(hit, isVsyn);
%         part = syn.pre(hit);
% 
%     else
%         isVsyn = [];
%         for v = 1:length(vgcCids)
%             isVsyn = [isVsyn; find(syn.pre == vgcCids(v))];
%         end
%         hit = intersect(hit, isVsyn);
%         part = syn.post(hit);
%     end
% 
%     if ~s(i).checkCids
%         s(i).synIds = hit;
%     else
%         isPart = part * 0;
%         for h = 1:length(hit)
%             if sum(s(i).cids == part(h))
%                 isPart(h) = 1;
%             end
%         end
%         s(i).synIds = hit(isPart>0);
%     end
% end

%% Check IDs

for i = 1:length(s)

    pre = syn.pre(s(i).synIds);
    post = syn.post(s(i).synIds);
    known = (pre > 0) & (post > 0);
    pre = pre(known);
    post = post(known);

    if ~isempty(pre)
        pairInd = sub2ind([max(pre) max(post)], pre,post);
        uPairInd = unique(pairInd);
        s(i).hPair = hist(pairInd,uPairInd);
    end
 
    if ~isempty(pre)
        uPairInd = unique(pre);
        s(i).hPre2All = hist(pre,uPairInd);
    end

     if ~isempty(post)
        uPairInd = unique(post);
        s(i).hPost2All = hist(post,uPairInd);
    end

    uPre = unique(pre);
    s(i).uPre = uPre;
    for p = 1:length(uPre)
        post1 = post(pre == uPre(p));
        uPost1 = unique(post1);
        hPost1 = histc(post1,uPost1);
        s(i).preP(p).cid = uPre(p);
        s(i).preP(p).uPost = uPost1;
        s(i).preP(p).hPost = hPost1;
    end

     uPost = unique(post);
     s(i).uPost = uPost;
    for p = 1:length(uPost)
        pre1 = pre(post == uPost(p));
        uPre1 = unique(pre1);
        hPre1 = histc(pre1,uPre1);
        s(i).postP(p).cid = uPost(p);
        s(i).postP(p).uPre = uPre1;
        s(i).postP(p).hPre = hPre1;
    end

end


%% Show hist

subN1 = round(sqrt(length(s)));
subN2 = ceil(length(s)/subN1);
hRange = [1:20];
clf
clear stdCon meanCon nCon resCon
for i = 1:length(s)
    subplot(subN1,subN2,i)
    useCount = s(i).hPre2All;
    %subplot(size(g2gN,1),size(g2gN,2),(g-1)*size(g2gN,1)+t)
    conCount = histc(useCount,hRange);
    meanCon(i) = mean(useCount);
    stdCon(i) = std(useCount)/sqrt(length(useCount));
    nCon(i) = length(useCount);
    bar(hRange,conCount/sum(conCount),'facecolor',[0 0 0],'edgealpha',0);
    title(sprintf('%s',s(i).name))
    %ylim([0 .2])
    drawnow
    resCon{i,1} = sprintf('%s %0.1f %0.1f %0.1f',s(i).name,nCon(i),meanCon(i),stdCon(i));
end


clf
useSG =5;
conPair = s(useSG).hPair;
conCount = histc(conPair,hRange);
bar(hRange,conCount/sum(conCount),'facecolor',[0 0 0],'edgealpha',0);
meanCon = mean(conPair)
n = length(conPair)
seCon = std(conPair)/sqrt(n)

drawnow


%% Show divergence

useSG =5;
clf

preP = s(useSG).preP;
postNum = zeros(length(preP),1);
synNum = postNum;
for p = 1:length(preP)
    synNum(p) = sum(preP(p).hPost);
    postNum(p) = length(preP(p).hPost);
end

meanPost = mean(postNum)
n = length(postNum)
sePost = std(postNum)/sqrt(n)


uSynNum = unique(synNum);
clear meanP maxP minP
for p = 1:length(uSynNum)
    
    synNums = postNum(synNum==uSynNum(p));
    meanP(p) = mean(synNums);
    se = std(synNums)/sqrt(length(synNums));
    highSE(p) = meanP(p) + se;
    lowSE(p) = meanP(p) - se;
    maxP(p) = max(synNums);
    minP(p) = min(synNums);

end

Y = [highSE fliplr(lowSE)];
X = [uSynNum(:)' fliplr(uSynNum(:)')];

postNumShow = postNum + (rand(length(postNum),1)-.5)/6;
scatter(synNum,postNumShow,'marker','o','markerfacealpha',.2,'markerfacecolor','k')
hold on
fill(X,Y,'r','facealpha',.2,'edgealpha',0)
plot(uSynNum,meanP,'linewidth',3,'color',[1 0 0])
hold off



%% Print figure
if 0
    %fDir = uigetdir;
    filename = [fDir '\bipConvergence_group5']
    set(gcf,'renderer','Painters')
    print('-depsc','-tiff','-r300', '-painters',[filename,'.eps'])

end










