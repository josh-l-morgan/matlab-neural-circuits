function[cons] = preTo(allEdges,targ)

isTarg = allEdges(:,2) == targ;

isCon = allEdges(isTarg,1);

uCon = unique(isCon);

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
    
    
    
    
    