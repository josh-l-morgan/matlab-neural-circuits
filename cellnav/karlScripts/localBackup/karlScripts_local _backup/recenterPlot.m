function recenterPlot(axs)
d=datacursormode();
inf=getCursorInfo(d);
[cursX,cursY,cursZ]=deal(inf.Position(1),inf.Position(2),inf.Position(3));
[curLimX,curLimY,curLimZ]=deal(axs.XLim,axs.YLim,axs.ZLim);
Xspan=max(abs(curLimX-cursX));
Yspan=max(abs(curLimY-cursY));
Zspan=max(abs(curLimZ-cursZ));
axs.XLim=[cursX-Xspan cursX+Xspan];
axs.YLim=[cursY-Yspan cursY+Yspan];
axs.ZLim=[cursZ-Zspan cursZ+Zspan];
end