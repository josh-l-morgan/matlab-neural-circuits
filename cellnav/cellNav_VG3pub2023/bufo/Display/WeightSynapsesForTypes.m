

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

%% set variables
useVgc = [2 3 4 5 13 14];
lc = 16;


%% colect data on cells

cellSynCount = zeros(length(s),length(useVgc));
clear synPos cellSynPos cellSource sClose synIdsA
countG = zeros(length(s),1);
for i = 1:length(s)
    synPos{i} = [];
    cellSource{i} = [];
    sClose{i} = [];
    synIdsA{i} = [];
end


%%  pick synapses

sNames = {s.name};
clear sg sNum


sg(1).names = {'bc3a'};
sg(1).col = [0 0 1];
sg(2).names = {'bc3b'};
sg(2).col = [0 .5 1];
sg(3).names = {'bc4'};
sg(3).col = [0 1 1];
sg(4).names = {'xbc'};
sg(4).col = [.5 1 0];
sg(5).names = {'bc5i'};
sg(5).col = [1 1 0];
sg(6).names = {'bc5o'};
sg(6).col = [1 .7 0];
sg(7).names = {'bc5t'};
sg(7).col = [1 .5 0];
sg(8).names = {'bc6'};
sg(8).col = [1 0 0];

sg(9).names = {'4on' '4i' '4ow'};
sg(9).col = [1 0 0];
sg(10).names = {'5ti'};
sg(10).col = [.6 .6 0];
sg(11).names = {'6sw'};
sg(11).col = [0 0 1];

doRGCs = 9:11;
doBips = 1:8;


%% pick s's for sgs
expNames = cat(1,{sg(:).names});
for g = 1:length(sg)
    nams = sg(g).names;
    for n = 1:length(nams)
        m1 = find(strcmp(sNames,sg(g).names{n}));
        if ~isempty(m1)
            sg(g).sT(n) = m1;
        else
            sg(g).sT(n) = [];
        end
    end
end



%% rand each cell

pol = zeros(length(sg),length(useVgc));
clear posV
for v = 1:length(useVgc)
    cid = useVgc(v);
    useSM(v) = 1;
    disp(sprintf('running sm %d, %d of %d',cid,v,length(useVgc)))

    targ = find(runCids==cid);
    sm = sms(targ).sm;%load([smDir fileName],'syn','syn2Skel','nep');

    syn = sm.syn; % get syn information
    synD = sm.syn2Skel.syn2SynDist;
    skelD = sm.syn2Skel.syn2SkelDist;
    D = sm.skel2skel.linDist; %use distance of every skeleton node to every other skeleton node
    nodeLengths = sm.nep.props.nodeLength;
    closest = sm.syn2Skel.closest;
    synIds = classifySynInCell(syn,s,cid);


    %allPos = cat(1,allPos,sm.arbor.subs);
    %%Get new node positions
    pos = sm.nep.pos;
    posV{v} = pos;

    %%Measure
    clear sgIds
    for g = 1:length(sg)
        sgIds{g} = cat(1,synIds{sg(g).sT});
        sgClose{g} = closest(sgIds{g});
        sNum(g,v) = length(sgIds{g});
    end
    clear minD meanD wD %% wD is influence of every group on every group
    for g1 = doRGCs % run for each rgc type
        for g2 = doBips
            if v== 1 %if first cell
                wGG{g2,g1} = [];
                posGG{g2,g1} = [];
            end
            d = D(sgClose{g2},sgClose{g1});
            W = exp(-d/lc);
            if sum(d(:))
                wG = sum(W,2); %sum synapses for rgc type
                sPos = syn.pos(sgIds{g2},:);
                wGG{g2,g1} = [wGG{g2,g1}; wG];
                posGG{g2,g1} = [posGG{g2,g1}; sPos];
            end
        end
    end


end


%% Show all
posA = [];
for v = 1:length(posV)
    posA = cat(1,posA,posV{v});
end

%%normalize
totGG = zeros(length(sg)); 
for g1 = doRGCs % run for each rgc type
        for g2 = doBips
        totGG(g2,g1) = sum (wGG{g2,g1});
        end
end
rgcTot = sum(totGG,1);

dotScale = 750
uDim = [1 3];
clf
c = 0;
for g1 = doRGCs % run for each rgc type
    c = c+1;
    subplot(3,1,c)
    hold on
%     scatter(posA(:,uDim(1)),posA(:,uDim(2)),3,'markerfacecolor',[.7 .7 .7],...
%         'MarkerEdgeColor','none')
    title(sg(g1).names)
    axis equal

    for g2 = doBips

        wG = wGG{g2,g1} / rgcTot(g1);
        wGSize = (wG * 100) * dotScale;
        sPos = posGG{g2,g1};
        scatter(sPos(:,uDim(1)),sPos(:,uDim(2)),wGSize,'markerfacecolor',sg(g2).col,...
            'markeredgecolor',sg(g2).col * .2,'markerfacealpha',.5)
    end

end


if 0
    epsName = 'Z:\Active\morganLab\PUBLICATIONS\VG3\Figure\Draft4\Pics3\BubbleInfluence5.eps';
    print(gcf, epsName, '-depsc2','-painters','-r300')
end


if 0
    clf
    hold on
    for g2 = doBips
        scatter([1 2 3],[g2 g2 g2],[0.01 .1 1] * dotScale,'markerfacecolor',sg(g2).col,...
            'markeredgecolor',sg(g2).col * .2,'markerfacealpha',.5)
    end
    xlim([0 4])
    ylim([min(doBips)-1 max(doBips)+1])
    epsName = 'Z:\Active\morganLab\PUBLICATIONS\VG3\Figure\Draft4\Pics3\BubbleInfluenceLegend5.eps';
    print(gcf, epsName, '-depsc2','-painters','-r300')
end

%%





