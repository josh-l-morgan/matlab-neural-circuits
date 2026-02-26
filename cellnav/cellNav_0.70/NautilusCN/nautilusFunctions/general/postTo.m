function[cons] = postTo(allEdges,targ)
%% function[cons] = postTo(allEdges,targ)


if isempty(allEdges)
    cons = [];
else

isTarg = allEdges(:,1) == targ;

isCon = allEdges(isTarg,2);

uCon = unique(isCon);
uCon = uCon(uCon>0);

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