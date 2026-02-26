function[picked] = pickWeighted(listIn)

sumList = cumsum(listIn.*(listIn>0))- 0.0000001;
picked =   sum(sumList<(rand*sumList(end)))+1;
if picked>length(listIn)
    picked = [];
end
