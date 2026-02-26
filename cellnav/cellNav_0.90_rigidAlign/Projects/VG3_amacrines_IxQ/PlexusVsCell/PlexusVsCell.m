
%%Run getSynapseNumberWithLocalExclusion.m to get cellCounts
runBipOrRgc = 1; %1 = RGCs
minDist = .01;



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
    load([SPN 'NOI.mat']);maPol
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


useVgc =[ 2 3 4 5 13 14];
for i = 1:length(sms)
    smCids(i) = sms(i).sm.cid;
end




bipCids = COI.bipCids;
clear bid2RGC posV
for v = 1:length(useVgc)
    cid = useVgc(v);
    useSM(v) = 1;
    disp(sprintf('running sm %d, %d of %d',cid,v,length(useVgc)))
    
    targ = find(runCids==cid);
    sm = sms(targ).sm;%load([smDir fileName],'syn','syn2Skel','nep');

    syn2syn = sm.syn2Skel.syn2SynDist;
    

    isToRGC = find(sm.syn.postClass == 1);

    pre = sm.syn.pre;

    syn2rgcSyn = syn2syn(:,isToRGC);
    posV{v} = sm.syn.pos(isToRGC,:);
    W =exp(-syn2rgcSyn/lc);

    bid2RGC{v} = zeros(length(bipCids),length(isToRGC));
    for b = 1:length(bipCids);
        isId = find(pre==bipCids(b));
        bid2RGC{v}(b,:) = sum(W(isId,:),1);
    end


end

b2r = [];
onV = [];
pos = [];
for v = 1:length(useVgc);
    b2rS = bid2RGC{v};
    b2r = cat(2,b2r,b2rS);
    onV = cat(2,onV,ones(1,size(b2rS,2))*v);
    pos = cat(1,pos,posV{v});
end

image(b2r)


%% analyze

b2rN = b2r ./ repmat(sum(b2r,1),[size(b2r,1) 1]);
image(b2rN*256)

bNum = size(b2rN,1);
rNum = size(b2rN,2);

bDist = zeros(rNum,rNum); % Distance in bipolar cell influence space
eDist = zeros(rNum,rNum); % Distance in euclidian distance
vDist = zeros(rNum,rNum,'logical'); % Same or different VG3
gDist = zeros(rNum,rNum,'logical');
for r1 = 1:rNum
    for r2 = 1:rNum
        dif = b2rN(:,r1) - b2rN(:,r2);
        bDist(r1,r2) = sqrt(sum(dif.^2));
        dif = pos(r1,:) - pos(r2,:);
        eDist(r1,r2) =  sqrt(sum(dif.^2));
        vDist(r1,r2) = onV(r1) ~= onV(r2);
        gDist(r1,r2) = r1~=r2;

    end
end

bDists = bDist(gDist);
eDists = eDist(gDist);
vDists = vDist(gDist);

sameCell = find(~vDists);
difCell = find(vDists);


for i = 1
    showSame = ceil(rand(300,1) * length(sameCell));
    showDif = ceil(rand(300,1) * length(difCell));
    clf, hold on
    scatter(eDists(difCell(showDif)),bDists(difCell(showDif)),'r')
    scatter(eDists(sameCell(showSame)),bDists(sameCell(showSame)),'g')
    drawnow
end


eRange = [0:1:150];

clear sameMeanB difMeanB sameSeB difSeB
for e = 1:length(eRange) -1
    isInRange = ((eDists>=eRange(e)) & (eDists<eRange(e+1)));
    vals = bDists(isInRange & ~vDists);
    sameMeanB(e) = mean(vals);
    sameSeB(e) = std(vals)/sqrt(length(vals));

    vals = bDists(isInRange & vDists);
    difMeanB(e) = mean(vals);
    difSeB(e) = std(vals)/sqrt(length(vals));
end


clf, hold on
eRange = eRange(1:end-1)

plot(eRange,sameMeanB,'g');
plot(eRange,sameMeanB-sameSeB,'g');
plot(eRange,sameMeanB+sameSeB,'g');

plot(eRange,difMeanB,'r');
plot(eRange,difMeanB-difSeB,'r');
plot(eRange,difMeanB+difSeB,'r');















