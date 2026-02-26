function[newSub] = downSampSub(rawSub,dSamp)

%% downSample object by dSamp


    sub = double(rawSub);
    if ~isempty(sub)
        sub(:,1) = round(sub(:,1)/dSamp(1));
        sub(:,2) = round(sub(:,2)/dSamp(2));
        sub(:,3) = round(sub(:,3)/dSamp(3));
        
        sub(sub<1) = 1; %%!!! Problem getting subs should never get zeros
        maxSub = max(sub,[],1);
        if ~isempty(sub)
            inds = sub2ind(maxSub,sub(:,1),sub(:,2),sub(:,3));
            uInds = unique(inds);
            if length(uInds)>1
               hInds = hist(inds,uInds);
               
            else
                hInds = length(inds)
            end
            [y x z] = ind2sub(maxSub,uInds);
            newSub = cat(2,uint16(y),uint16(x),uint16(z));
            %dsObj(i).n = uint16(hInds);
        end
    else %if sub is empty
        newSub = rawSub;
    end
    