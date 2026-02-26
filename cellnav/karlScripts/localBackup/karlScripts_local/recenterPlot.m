function recenterPlot(axs)
% if ~exist('labLoc','var')
%     labLoc='NON'
% end
d=datacursormode();
inf=getCursorInfo(d);
[cursX,cursY,cursZ]=deal(inf.Position(1),inf.Position(2),inf.Position(3));
[curLimX,curLimY,curLimZ]=deal(axs.XLim,axs.YLim,axs.ZLim);
Xspan=max(abs(curLimX-cursX));
Yspan=max(abs(curLimY-cursY));
Zspan=max(abs(curLimZ-cursZ));
axs.XLim=[cursX-(Xspan/2) cursX+(Xspan/2)];
axs.YLim=[cursY-(Yspan/2) cursY+(Yspan/2)];
axs.ZLim=[cursZ-(Zspan/2) cursZ+(Zspan/2)];
end