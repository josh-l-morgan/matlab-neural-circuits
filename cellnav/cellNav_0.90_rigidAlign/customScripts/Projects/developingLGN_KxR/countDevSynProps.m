
pFig = figure;
pAx = gca(pFig);
pAx.NextPlot = "add";

sp = parseDevSyn();
isCid = tis.syn.post == 3097;

rgc = sp.tp.rgc;
lin = sp.tp.lin;

sp.st.preClass.char;  %%Tag order for syn properties
preSynVal = sp.ps.preSynVal;


rgcSheath = rgc & sp.ps.isSheath;
rgcNoSheath = rgc & ~sp.ps.isSheath;

numRgcSheath = sum(rgcSheath & isCid);
numRgcNoSheath = sum(rgcNoSheath & isCid);

clear isSum
isType = sp.ps.isSheath;
for i = 1:4
    isBut = rgc & isCid & (preSynVal(:,1)==i);
    isSum(1,i) = sum(isBut & isType);
    isSum(2,i) = sum(isBut & ~isType);
end

plot(pAx,isSum(1,:),'r');
plot(pAx,isSum(2,:),'b');











