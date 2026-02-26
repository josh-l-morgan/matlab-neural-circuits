function[cons] = preTo(allEdges,targ)

if isempty(allEdges)
    cons = [];
else
isTarg = allEdges(:,2) == targ;

isCon = allEdges(isTarg,1);

uCon = unique(isCon);

if length(uCon)>1;
    hCon = histc(isCon,uCon);
else
    hCon = length(isCon);
end
   
if sum(hCon)
    cons  = [uCon hCon];
else
    cons = [];
end
end
    
    
    
    
    