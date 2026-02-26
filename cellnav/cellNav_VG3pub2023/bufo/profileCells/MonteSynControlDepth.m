
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
    
clear useVgc
for i = 1:length(sms)
    useVgc(i) = sms(i).sm.cid;
end



%% colect data on cells

cellSynCount = zeros(length(s),length(useVgc));
clear synPos cellSynPos cellSource sClose
countG = zeros(length(s),1);
for i = 1:length(s)
    synPos{i} = [];
    cellSource{i} = [];
    sClose{i} = [];
end
allPos = [];
for v = 1:length(useVgc)

    cid = useVgc(v);
    useSM(v) = 1;
    disp(sprintf('loading sm %d, %d of %d',cid,v,length(useVgc)))
    sm = sms(v).sm;%load([smDir fileName],'syn','syn2Skel','nep');

    syn = sm.syn; % get syn information
    synD = sm.syn2Skel.syn2SynDist;
    skelD = sm.syn2Skel.syn2SkelDist;
    nodeLengths = sm.nep.props.nodeLength;  
    closest = sm.syn2Skel.closest;


    allPos = cat(1,allPos,sm.arbor.subs);

    %%Create synapse groups filled with relevant IDs using COI and subtypes


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
            synIds = hit;
        else
            isPart = part * 0;
            for h = 1:length(hit)
                if sum(s(i).cids == part(h))
                    isPart(h) = 1;
                end
            end
            synIds = hit(isPart>0);
        end
        cellSynCount(i,v) = length(synIds);

        cellSynPos{i,v} = syn.pos(synIds,:);

        synPos{i} = cat(1,synPos{i},syn.pos(synIds,:));
        cellSource{i} = cat(1,cellSource{i},synIds *0+v);
        sClose{i} = cat(1,sClose{i}, closest(synIds));
    end


end

allSynCount = sum(cellSynCount,2);


%%  pick synapses

sNames = {s.name};
clear sg
sg(1).names = {'off'};
sg(1).col = [0 1 0];
sg(2).names = {'bc6'};
sg(2).col = [1 0 0];
sg(3).names = {'4ow'};
sg(3).col = [0 0 1];
sg(4).names = {'51'};
sg(4).col = [.6 .6 0];


%% pick s's for sgs
expNames = cat(1,{sg(:).names});
for g = 1:length(sg)
    nams = sg(g).names;
    sg(g).pos = [];
    sg(g).v = [];
    for n = 1:length(nams)
        m1 = find(strcmp(sNames,sg(g).names{n}));
        if ~isempty(m1)
            sg(g).sT(n) = m1;
            sg(g).pos = cat(1,sg(g).pos,synPos{m1});
            sg(g).v = cat(1,sg(g).v,cellSource{m1});
        else
            sg(g).sT(n) = [];
        end
    end

end


%% transform synapses
clear synDepth
for i = 1:length(sg)
    pos = sg(i).pos;
    [zGCL,zINL,sg(i).depth] = getIPLdepth(pos(:,3),pos(:,1),pos(:,2),[],[]);
end


[zGCL,zINL,nDepth] = getIPLdepth(allPos(:,3),allPos(:,1),allPos(:,2),[],[]);



%% Find depths
%binDepths(sg,nDepth)

%% rand each cell
zJit = .5;
reps = 10;

for v = 1:length(useVgc)
    cid = useVgc(v);
    useSM(v) = 1;
    disp(sprintf('loading sm %d, %d of %d',cid,v,length(useVgc)))
    sm = sms(v).sm;%load([smDir fileName],'syn','syn2Skel','nep');

    syn = sm.syn; % get syn information
    synD = sm.syn2Skel.syn2SynDist;
    skelD = sm.syn2Skel.syn2SkelDist;
    nodeLengths = sm.nep.props.nodeLength;
    closest = sm.syn2Skel.closest;


    %allPos = cat(1,allPos,sm.arbor.subs);
    %%Get new node positions
    pos = sm.nep.pos;
    [zGCL,zINL,nDepth] = getIPLdepth(pos(:,3),pos(:,1),pos(:,2),[],[]);
    iplDepth = mean(zINL)-mean(zGCL);
    nDepth = nDepth * iplDepth;

    zDif = abs(nDepth-nDepth');
    zMask = zDif<zJit;

    clear isV
    for i = 1:length(sg)
        isV{i} = find(sg(i).v ==v);
    end
    


    for r = 1:reps

        %make random pointers from nodes to node with similar z
        while 1
            zRand = rand(size(zMask)) .* zMask;
            zMax = repmat(max(zRand,[],1),[size(zRand,1) 1]);
            [rClose x] = find(zRand==zMax);
            if length(rClose) == size(zMax,1)
                break
            end
        end


        for i = 1:length(sg)

            oldClose = sClose{i}(isV{i});
            newClose = rClose(oldClose);
            if 1
                clf

                isBip = find(syn.preClass == 7);

                scatter3(pos(:,1),pos(:,2),pos(:,3),'k','.');
                hold on
                %scatter3(pos(closest(isBip),1),pos(closest(isBip),2),pos(closest(isBip),3),'m','o','filled');

                scatter3(pos(oldClose,1),pos(oldClose,2),pos(oldClose,3),'r','filled');
                scatter3(pos(newClose,1),pos(newClose,2),pos(newClose,3),'g','filled');
                hold off
                axis equal
                set(gca,'clipping', 'off')
                view([0 0])
            end
            pause
        end

    end




end
%{
function binDepths(sg,nDepth)

%% binDepths
clf
bin = .04;
dRange = -.1:.001:1.1;

maxCount = 0;

for i = 1:length(sg)

    dCount = dRange*0;
    for b = 1:length(dRange)
        d = dRange(b);
        dCount(b) = sum((sg(i).depth>(d-bin/2)) & (sg(i).depth<=(d+bin)));
    end
    sg(i).dCount = dCount;
    maxCount = max(maxCount,max(dCount));
end

clf
subplot(2,1,1)
hold on
for i = 1:length(sg)
    plot(dRange,sg(i).dCount/maxCount,'color',sg(i).col)
end

subplot(2,1,2)
hold on
nCount = dRange*0;
bin2 = .1;
for b = 1:length(dRange)
    d = dRange(b);
    nCount(b) = sum((nDepth>(d-bin2/2)) & (nDepth<=(d+bin2)));
end
nCount = nCount/max(nCount(:));
plot(dRange,nCount,'color',[.3 .3 .3])


ei = sg(1).dCount./sg(2).dCount;
plot(dRange,ei,'color',[0 0 1])



iDens = (sg(2).dCount/maxCount)./(nCount/max(nCount(:)));
plot(dRange,iDens,'color',[1 0 0])
% 
% 
% eDens = (sg(1).dCount/maxCount)./(nCount/max(nCount(:)));
% plot(dRange,eDens,'color',[0 1 0])

end

%}






