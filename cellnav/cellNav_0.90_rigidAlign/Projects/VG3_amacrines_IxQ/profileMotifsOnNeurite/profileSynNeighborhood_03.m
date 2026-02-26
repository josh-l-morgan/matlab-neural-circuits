%Generate violyn scatter plots of nearby synapses. Grouped and colored by type.
 
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


%% set variables
distBin = 3;
bRange = [0 : .1 : 20];
ignoreSelf = 0; %ignore synapses between the same sets of cell IDs, 2 = only self




%% colect data on cells
g2g = zeros(length(s),length(s),length(bRange));
nearLength = g2g;
countG = zeros(length(s),1);
for v = 1:length(useVgc)

    cid = useVgc(v);
    useSM(v) = 1;
    disp(sprintf('loading sm %d, %d of %d',cid,v,length(useVgc)))
    sm = sms(v).sm;%load([smDir fileName],'syn','syn2Skel','nep');

    syn = sm.syn; % get syn information
    synD = sm.syn2Skel.syn2SynDist;
    skelD = sm.syn2Skel.syn2SkelDist;
    nodeLengths = sm.nep.props.nodeLength;


    synD(sub2ind(size(synD),1:size(synD),1:size(synD))) = inf; %set distance to self as infinite

    if ignoreSelf %% ignore synapses between the same pair of cells
        for sY = 1:size(synD,1)
            pre1 = syn.pre(sY);
            post1 = syn.post(sY);
            hit = ((syn.pre==pre1) & ( syn.post == post1))>0;
            if ignoreSelf == 1 % remove same cell pairings
                synD(sY,hit) = inf;
            elseif ignoreSelf == 2 % remove everything but same cell pairings
                synD(sY,~hit) = inf;
            end
        end
    end


    %% sort types
    synType = zeros(length(syn.preClass),1);
    synType( (syn.preClass == 7) | (syn.synType == 2) ) = 1; % bip
    synType( (syn.synType == 1) & (syn.post == cid) ) = 2; % amc in
    synType( (syn.preClass == 0) & (syn.post == cid) ) = 2; % amc in
    synType( (syn.pre == cid) & (syn.postClass == 1) ) = 3; % rgc
    synType( (syn.pre == cid) & (syn.postClass == 8) ) = 4; % amc out
    syn.post(synType == 0);

    %% Create synapse groups filled with relevant IDs using COI and subtypes


    %%Collect synapses
    for i = 1:length(s)
        if ~isempty(s(i).synType)
            hit = find(syn.synType == s(i).synType);
        else
            hit = find(syn.synType<inf) ;
        end

        if s(i).input
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
        end
    end


    %% Bin presence of different synapses and lengths

    for g = 1:length(s);
        disp(sprintf('checking group %s, %d of %d',s(g).name,g,length(s)))
        sG = s(g).synIds;
        countG(g) = countG(g) + length(sG);
        synDG = synD(sG,:);
        skelDG = skelD(sG,:);
        for hT = 1:length(s) %check against all types
            synDG2 = synDG(:,s(hT).synIds);
            for b = 1:(length(bRange))
                [ y x] = find((synDG2>=(bRange(b)-distBin/2)) & (synDG2<(bRange(b)+distBin/2)));
                g2g(g,hT,b) =  g2g(g,hT,b) + length(x);
            end
        end
        for b = 1:(length(bRange)) %get length for density
            [ y x] = find((skelDG>=(bRange(b)-distBin/2)) & (skelDG<(bRange(b)+distBin/2)));
            nearLength(g,:,b) = nearLength(g,:,b) + sum(nodeLengths(x));
        end

    end

end

%% Scale by number of gs found
g2g = g2g ./ repmat(countG,[1 size(g2g,2) size(g2g,3)]);
nearLength = nearLength ./ repmat(countG,[1 size(g2g,2)]);

g2gN = g2g ./ nearLength;


%% Display group to group
if 0
    clf
    hold on
    tCol = [1 0 0; 0 1 0 ; 0 0 1; 1 1 0; 1 0 1; 0 1 1]*.7;
    sumT = squeeze(sum(g2g,2));
    maxT = max(sumT(:));

    subN1 = round(sqrt(length(s)));
    subN2 = ceil(length(s)/subN1);
    for g = 1:size(g2g,1)
        clf
        set(gcf,'Name',s(g).name)
        for t = 1:size(g2g,2)
            subplot(subN1,subN2,t)
            %subplot(size(g2gN,1),size(g2gN,2),(g-1)*size(g2gN,1)+t)
            plot(bRange,squeeze(g2gN(g,t,:)))
            title(sprintf('%s',s(t).name))
            ylim([0 .25])

        end
        drawnow
        pause
    end

    pause(.1)

end



%% Display select groups
clf
clear nh
sNames = {s(:).name}'
if 0
    nh(1).ref = '4 4i 4on 4ow';
    nh(1).check = {'bc3a' 'bc3b' 'bc4' 'bc5i' 'bc5t' 'bc6' 'non-rib'};

    nh(2).ref = '4ow';
    nh(2).check = {'bc3a' 'bc3b' 'bc4' 'bc5i' 'bc5t' 'bc6' 'non-rib'};

    nh(2).ref = '4i';
    nh(2).check = {'bc3a' 'bc3b' 'bc4' 'bc5i' 'bc5t' 'bc6' 'non-rib'};

    nh(3).ref = '51';
    nh(3).check = {'bc3a' 'bc3b' 'bc4' 'bc5i' 'bc5t' 'bc6' 'non-rib'};

    nh(4).ref = '37 37c 37d 37r 37v';
    nh(4).check = {'bc3a' 'bc3b' 'bc4' 'bc5i' 'bc5t' 'bc6' 'non-rib'};

    nh(5).ref = '63';
    nh(5).check = {'bc3a' 'bc3b' 'bc4' 'bc5i' 'bc5t' 'bc6' 'non-rib'};

    nh(6).ref = 'RGCs';
    nh(6).check = {'bc3a' 'bc3b' 'bc4' 'bc5i' 'bc5t' 'bc6' 'non-rib'};
else


    nh(1).ref = 'rib' ;
    nh(1).check = {'rib' 'RGCs' 'non-rib'};

    nh(2).ref = 'RGCs';
    nh(2).check = {'rib' 'RGCs' 'non-rib'};

    nh(3).ref = 'non-rib' ;
    nh(3).check = {'rib' 'RGCs' 'non-rib'};




end


checkCol = [0 0 1; 0 .5 .5; 0 1 1; 1 0 0; 1 .5 0; .5 .5 0; 0 0 0];
clf
for i = 1:length(nh)
    subplot(ceil(length(nh)/2),2,i)
    hold on
    m1 = find(strcmp(sNames,nh(i).ref))
    check = nh(i).check;
    titleStr = sprintf('%s, n = %d',nh(i).ref,countG(m1))
    title(titleStr)
    for c = 1:length(check)
        m2 = find(strcmp(sNames,check{c}))
        m12 = squeeze(g2gN(m1,m2,:));
        if sum(checkCol(c,:))
            plot(bRange,m12,'color',checkCol(c,:));
        else
            plot(bRange,m12/4,'color',checkCol(c,:));
        end
    end
    hold off
    ylim([0 .25])
    if i == length(nh)
        legend(nh(i).check)
    end
end

return

%% Scatter syn
clf
hold on

for t = 1:size(g2gN,2)
    scatStore(t).X = [];
    scatStore(t).Y = [];
end

g2gN2 = g2gN * 100;

spacer = .6/max(g2gN2(:));
for g =  1:size(g2gN2,1)

    for b = 1:min(10,size(g2gN2,3))
        tVec = g2gN2(g,:,b);
        sumVec = sum(tVec);
        for t = 1:length(tVec)
            if tVec(t)
                n = tVec(t);
                startX = (n*spacer)/2 * -1;
                shiftX = [startX:spacer:startX+(n * spacer)] + g;
                shiftY = (shiftX*0+bRange(b)) + ((t/length(tVec))*distBin*.6);
                %scatter(shiftX,bRange(b),10,tCol(t,:))
                scatStore(t).X = cat(2,scatStore(t).X,shiftX);
                scatStore(t).Y = cat(2,scatStore(t).Y,shiftY);

            end
        end

    end
end
for t = 1:size(g2gN2,2)

    scatter(scatStore(t).X,scatStore(t).Y,30,tCol(t,:),'filled')
    hold on

end

return

%% RGC to BPC
clf
plot(squeeze(g2gN(3,1,:)),'g')
hold on
plot(squeeze(g2gN(3,3,:)),'b')
hold off
ylim([0 .15])


%% Print figure
if 0
    %fDir = uigetdir;
    filename = [fDir '\rgcNeighborhoods']
    set(gcf,'renderer','Painters')
    print('-depsc','-tiff','-r300', '-painters',[filename,'.eps'])

end







