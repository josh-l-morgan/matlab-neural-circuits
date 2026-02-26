function[nep] = smoothRad(nep,spread)


rawRad = nep.nodeRad;
%
% edgeRef = nep.edges(:);
% uNode = unique(edgeRef);
% hNode = hist(edgeRef,nep.nodes);

sumNodeVal = zeros(length(nep.nodes),1);
countNodeVal =zeros(length(nep.nodes),1);
saveVal = zeros(length(nep.nodes),length(nep.bones));

for b = 1:length(nep.bones)
    edges = nep.bones(b).edges;
    eRun = cat(1,edges(1,1),edges(:,2));
    lengths = sqrt((nep.pos(edges(:,1),1)-nep.pos(edges(:,2),1)).^2 + ...
        (nep.pos(edges(:,1),2)-nep.pos(edges(:,2),2)).^2 + ...
        (nep.pos(edges(:,1),3)-nep.pos(edges(:,2),3)).^2);
    vals = rawRad(eRun);
    meanVal = zeros(length(eRun),1);
    for n = 1:length(eRun);
        
        L = 0;
        bRad = [];
        for e = 1:100
            if ((n-e)>0)
                L = L + lengths(n - e);
                if L< spread
                    bRad = cat(1,bRad,vals(n-e));
                else
                    break
                end
            else
                break
            end
        end
        
        L = 0;
        fRad = vals(n);
        for e = 1:100
            if ((n+e)<length(eRun))
                L = L + lengths(n + e-1);
                if L< spread
                    fRad = cat(1,fRad,vals(n+e));
                else
                    break
                end
            else
                break
            end
        end
        
        meanVal(n) = mean([bRad;fRad]);
    end

    nep.bones(b).meanRads = meanVal;
    
    
    sumNodeVal(eRun)  = sumNodeVal(eRun) + meanVal;
    countNodeVal(eRun) = countNodeVal(eRun) + 1;
    saveVal(eRun,b) = meanVal;
end

nep.meanNodeRad = max(saveVal,[],2);
    
    
    
    




