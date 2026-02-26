
%%Display synapse distributions of targeted TCs

runCids = COI.targetTCs;
%runCids = 3204;
numRun = length(runCids);

smCids = zeros(length(sms),1);
for i = 1:length(sms)
    try
        smCids(i) = sms(i).sm.cid;
    end
end

maxDist = 0;

synPreClasses = [0 1 2 3];
synClassColor = [.3 .3 .3; 0 1 0; 0 0 1; 1 0 1];
synClassAlpha = [.3 .3 .3 .8];
synClassSize = [50 50 50 100];

if 1
    fig = figure;
    ax(1) = subplot(2,1,1,'parent',fig);
    ax(1).NextPlot = 'add';
    ax(2) = subplot(2,1,2,'parent',fig);
    ax(2).NextPlot = 'add';
end

for i = 1:length(runCids)

    %%Fetch data
    smTarg = find(smCids==runCids(i),1);
    if ~isempty(sms(smTarg))
        sm = sms(smTarg).sm;
        nep = sm.nep; %light weight representation of skeleton and synapses
        syn = sm.syn; %synapse information

        %%Get distances to seeds
        seedNode = nep.seedNode; %Should be center of cell body
        D = sm.skel2skel.linDist; %matrix of all node to node distances;
        preClass = syn.preClass; %cell class of presynaptic partner
        closest = sm.syn2Skel.closest;

        node2SeedDistance = D(:,seedNode); %Find distance of all nodes to seed node
        syn2SeedDistance = node2SeedDistance(closest);

        %%Make topoGraph
        swc = nep.swcS;
        swc2SeedDistance = node2SeedDistance(swc.swc2arborID);
        maxDist = max(maxDist,max(swc2SeedDistance));
        topo = swc2topo(swc.pred,swc2SeedDistance);

        %%Plot branch graph
        cla(ax(2))
        ax(2).NextPlot = 'add';
        plot(ax(2),topo.graph.bpX,topo.graph.bpY,'color',[.9 .9 .9])
        plot(ax(2),topo.graph.branchX,topo.graph.branchY,'k');

        %%Add synapses
        closestSwc = swc.arbor2swcID(closest);
        syn2branch = topo.branch.IDofAllNodes(closestSwc);
        for c = 1:length(synPreClasses)
            showSyn = find(preClass==synPreClasses(c));
            scatter(ax(2), topo.graph.branchX(1,syn2branch(showSyn)), syn2SeedDistance(showSyn),'markerfacecolor',...
                synClassColor(c,:),'markerfacealpha',synClassAlpha(c),'MarkerEdgeColor',[0 0 0]);
        end

        %% swarm polot synapses
        lookupSynClassID(synPreClasses+1) = 1:length(synPreClasses);
        synClassDisplayID = lookupSynClassID(preClass+1);
        synCol = synClassColor(synClassDisplayID,:);
        synAlph = synClassAlpha(synClassDisplayID);
        synSize = synClassSize(synClassDisplayID);
        showSyn = 1:length(preClass);

        xJit = 10;
        sc = swarmchart(ax(2), showSyn * 0 + min(topo.graph.branchX(:))-xJit, syn2SeedDistance(showSyn),50,...
            'markerfacecolor','flat','markerfacealpha',.4,'MarkerEdgeColor',[0 0 0],'XJitter','density',...
            'XJitterWidth',xJit,'markerfacealpha',.3);
        sc.CData = synCol;
        sc.AlphaData = synAlph;
        sc.SizeData = synSize;
        title(ax(2),sprintf('topoGraph for cid%d',runCids(i)))



        sc2 = swarmchart(ax(1), showSyn * 0 + i * xJit*2, syn2SeedDistance(showSyn),50,...
            'markerfacecolor','flat','markerfacealpha',.4,'MarkerEdgeColor',[0 0 0],'XJitter','density',...
            'XJitterWidth',xJit,'markerfacealpha',.3);
        sc2.CData = synCol;
        sc2.AlphaData = synAlph;
        sc2.SizeData = synSize;

       


        pause(1)



    end

end

%%Set equal Ylim for all axis
ax(1).YLim = [0 maxDist+maxDist*.03];
for i = 1:length(runCids)
text(ax(1),i*xJit*2-xJit/2,maxDist + maxDist*.01,sprintf('%d',runCids(i)))
end






