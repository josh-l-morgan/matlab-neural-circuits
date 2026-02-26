


bX = 1 * 10^4;
bY = 2 * 10^4;
bZ = 1 * 10^4;

xCut = 3 * 10^3;
zCut = 2 * 10^1;
pCut = 2 * 10^1;


yStrips = bX / xCut
zStrips = bZ / zCut
pillars = xCut/pCut

stripTime = bY /500

stripNum = yStrips * zStrips * pillars

seconds = stripNum * stripTime

years = seconds/ 60 /60/24/365