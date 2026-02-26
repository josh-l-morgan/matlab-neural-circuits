%% diagnostic
f1=figure(); histogram(curSynDists,'BinEdges',[0:1:10]); title('distance');
xlim([0 10]); ylim([0 10]);

f2=figure(); histogram(curSynInfluence,'BinEdges',[0:.05:10]); title('influence');
xlim([.3 1]); ylim([0 10]);