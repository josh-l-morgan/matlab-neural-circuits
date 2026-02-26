

%% load nep

nepName = 'nepTest_128.mat';
nepDir = 'D:\LGNs1\Analysis\LIN\nepDir\';
load([nepDir nepName])
nep = nepCSS;

springRes = 'D:\LGNs1\Analysis\springDat\results\';
%load([springRes 'firstNepCSS.mat'])
load([springRes 'fixed125_2.mat'])

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

springDat = nep2SpringParameters_02(nep);


%% Run springs
%movDir  = [springDir 'seedMovie_?\'];
%if ~exist(movDir,'dir'), mkdir(movDir), end

for rerun = 1: 1
    allResults{rerun} = runSkelSprings(springDat);

end

%%  print eps
if 0
    springDir = 'D:\LGNs1\Analysis\springDat\skelSpring\';
    if ~exist(springDir,'dir'), mkdir(springDir), end
    
    set(gcf,'PaperUnits','points','PaperPosition',[1 1 700 700])
    set(gcf, 'InvertHardCopy', 'off');
    tag = 'net125_whitefull2';
    
    epsName = sprintf('%sspringRun_%s.eps',springDir,tag);
    print(gcf, epsName, '-depsc2','-painters','-r300')
end



%% Save result
if 0
    
    if ~exist(springRes,'dir'), mkdir(springRes), end
    
    result = allResults{end};
    save([springRes 'fixed125_2.mat'],'result')
    
    
end






