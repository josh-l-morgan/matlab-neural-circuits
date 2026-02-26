

%% load nep
load('MPN.mat')
nepName = 'nepTest_128b.mat';
nepDir = [WPN '\LIN\nepDir\'];
load([nepDir nepName])
nep = nepCSS;


%% Select from Nep

fNames = forceNodes();
fNodes = find(ismember(nep.nodeName,fNames));


%%Select nodes
startNodes = nep.nodes(nep.nodeType == 3);
linkedNodes = nearNodes(nep.edges,startNodes,2);
useNodes = unique([startNodes(:); linkedNodes(:); fNodes(:)]);
nep = trimNep(nep,useNodes);

%%Reshape Nep





%% grap data

nodeIDs = nep.nodes;
allEdges = nep.edges;
nodeNum = length(nodeIDs);
nodeType = nep.nodeType;
seedList = [ 108 201 907 903];

%% Set color

nep = colorNep(nep);


%% Create springDat (need nodeCol, nodeIDs,

nep.allWeights = nep.edges;
nep.seedList = seedList;

springDat = nep2SpringParameters_near(nep);


% Run springs
% movDir  = [springDir 'pull125_c\'];
% if ~exist(movDir,'dir'), mkdir(movDir), end

for rerun = 1: 1
    allResults{rerun} = runSkelSprings(springDat);

end

%%  print eps
if 0
    if ~exist(springDir,'dir'), mkdir(springDir), end
    
    set(gcf,'PaperUnits','points','PaperPosition',[1 1 700 700])
    set(gcf, 'InvertHardCopy', 'off');
    tag = 'net125_tempPart4';
    
    epsName = sprintf('%sspringRun_%s.eps',springDir,tag);
    print(gcf, epsName, '-depsc2','-painters','-r300')
end



%% Save result
if 0
    
    if ~exist(springRes,'dir'), mkdir(springRes), end
    
    result = allResults{end};
    save([springRes 'fixed125_2.mat'],'result')
    
    
end






