function plotsOut=plotSyns(synList,curTis,options)
%% plotSyns
% plot some synapses onto the compareMorph figures
% plotSyns(synList,curTis,options)
%global curTis;
arguments
       synList
       curTis
       options.depth = 0
       options.color = rand(1,3)
       options.size = 50
end
color=options.color;
depth=options.depth;
size=options.size;

if depth==1
    [~,~,synD]=getIPLdepth(curTis.syn.pos(synList,3),curTis.syn.pos(synList,1),curTis.syn.pos(synList,2),[],[]);
    synD=synD*40;
else
    synD=curTis.syn.pos(synList,3);
end

plotsOut=scatter3(synD,curTis.syn.pos(synList,1),curTis.syn.pos(synList,2), ...
    size, color);
end

