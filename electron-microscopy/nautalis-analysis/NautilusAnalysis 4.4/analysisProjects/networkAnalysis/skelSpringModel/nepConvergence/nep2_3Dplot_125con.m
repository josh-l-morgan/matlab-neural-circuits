

%% load nep
load('MPN.mat')
nepName = 'nepTest_128b.mat';
nepDir = [WPN '\LIN\nepDir\'];
load([nepDir nepName])
nep = nepCSS;

springRes = [WPN '\springDat\results\'];
%load([springRes 'firstNepCSS.mat'])
load([springRes 'fixed125_2.mat'])
    springDir = [WPN 'springDat\skelSpring\'];

%%%%%%%fake Nep
% 
% clear nep
% nep.nodes = 1:2;
% nep.edges = [1 2];%; 2 3; 3 4; 5 6;5 7; 5 8; 8 9; 9 10];
% nep.nodeType = ones(length(nep.nodes),1);

% 
% nep.nodes = 1:10;
% nep.edges = [1 2; 2 3; 3 4; 5 6;5 7; 5 8; 8 9; 9 10; 2 8];
% nep.nodeType = [1 1 1 2 2 2 3 3 3 3];
% nep.edgeType = [1 1 1 1 1 1 1 1 1];


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

nep = colorNep_con(nep);


%% Create springDat (need nodeCol, nodeIDs,

nep.allWeights = nep.edges;
nep.seedList = seedList;

springDat = nep2SpringParameters_near(nep);


plot3nep(nep)

% Run springs
% movDir  = [springDir 'pull125_c\'];
% if ~exist(movDir,'dir'), mkdir(movDir), end
% 
% for rerun = 1: 1
%     allResults{rerun} = runSkelSprings(springDat);
% 
% end






%%  print eps
if 0
    if ~exist(springDir,'dir'), mkdir(springDir), end
    
    set(gcf,'PaperUnits','points','PaperPosition',[1 1 700 700])
    set(gcf, 'InvertHardCopy', 'off');
    tag = 'numSynSpread20a';
    
    epsName = sprintf('%sspringRun_%s.eps',springDir,tag);
    print(gcf, epsName, '-depsc2','-painters','-r300')
end



%% Save result
if 0
    
    if ~exist(springRes,'dir'), mkdir(springRes), end
    
    result = allResults{end};
    save([springRes 'fixed125_2.mat'],'result')
    
    
end






