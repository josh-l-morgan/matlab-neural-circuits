clear all
%% Load data
load('MPN.mat')
%MPN = 'D:\LGNs1\Export\export_KV_LIN_morph_2019+3+7K\';

if ~exist('MPN','var')
    MPN = GetMyDir
end

synDir = [MPN 'synPos3\'];
if ~(exist(synDir,'dir')),mkdir(synDir); end

if ~exist('obI','var') | ~exist('dsObj','var')
    disp('loading')
    load([MPN 'obI.mat'])
    load([MPN 'dsObj.mat'])
end

cellList = {4};
subCell = names2Subs(obI,dsObj,cellList);
sub = subCell{1};
    
 clf
 renderCon(sub,[],[1 1 0],.5)
    
    
%% Translate

X = 1053; %<
Y = 1497; %>
Z = 403;  %><

anc2sub = (obI.em.res / 1000)./ obI.em.dsRes 

allPos2 = [X Y Z] ./ anc2sub;

sprintf('%.0f %.0f %.0f',allPos2(end,2),allPos2(end,1),allPos2(end,3))


