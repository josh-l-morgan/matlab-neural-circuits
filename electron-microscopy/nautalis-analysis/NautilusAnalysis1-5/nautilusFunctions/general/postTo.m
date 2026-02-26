function[cons] = postTo(allEdges,targ)

isTarg = allEdges(:,1) == targ;

isCon = allEdges(isTarg,2);

uCon = unique(isCon);
uCon = uCon(uCon>0);

if length(uCon)>1;
    hCon = hist(isCon,uCon);
else
    hCon = length(isCon);
end
    
    
if sum(hCon)
    cons  = [uCon hCon'];
else
    cons = [];
end
    