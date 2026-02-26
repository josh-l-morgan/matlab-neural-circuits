function[syn sProps] = getPrePost(allS,allO);


%%
labSyn = bwlabeln(allS,6);
sProps = regionprops(labSyn,allO,'PixelIdxList','PixelValues');
clear labSyn

%%
binOb = double([0:max(allO(:))]);
%allSyn = zeros(length(sProps),2);
numSyn = 0;
clear syn
for i = 1:length(sProps)
    vals = double(sProps(i).PixelValues);
    histVals = hist(vals,binOb);
    %histVals = histVals(2:end);
    hits = find(histVals>0);
    
    if length(hits)>=2
        
        hitVal = histVals(hits);
        [sortVal hIdx] = sort(hitVal,'descend');
        
        sProps(i).pre = hits(hIdx(1));
        sProps(i).post = hits(hIdx(2));
        numSyn = numSyn+1;
        syn(numSyn,1:2) = [hits(hIdx(1)) hits(hIdx(2))];
        synId(numSyn) = i;
    elseif length(hits) <2
        sProps(i).pre = hits;
        sProps(i).post = hits;
    end
    
end
