function[dsObj] = downSampObj(MPN, dSamp)

%% downSample object by 8

%MPN = 'D:\LGNs1\Segmentation\VAST\S8\joshm\export+14+04+27_mat\'
disp('loading vastSubs')
load([MPN 'vastSubs.mat'])
if ~exist('dSamp','var')
    dSamp = [1 1 1];
end

%% down sample

for i = 1:length(vastSubs)
    
    disp(sprintf('%d of %d',i,length(vastSubs)));
   
    dSampU = dSamp(min(size(dSamp,1),i),:);
    
    sub = double(vastSubs{i});
    if ~isempty(sub)
        sub(:,1) = double(vastSubs{i}(:,1))/dSampU(1);
        sub(:,2) = double(vastSubs{i}(:,2))/dSampU(2);
        sub(:,3) = double(vastSubs{i}(:,3))/dSampU(3);
        sub = round(sub);
        sub(sub<1) = 1; %%!!! Problem getting subs should never get zeros
        maxSub = max(sub,[],1);
        if ~isempty(sub)
            inds = sub2ind(maxSub,sub(:,1),sub(:,2),sub(:,3));
            uInds = unique(inds);
            if length(uInds)>1
               hInds = hist(inds,uInds);
            else
                hInds = length(inds);
            end
            [y x z] = ind2sub(maxSub,uInds);
            dsObj(i).subs = cat(2,uint16(y),uint16(x),uint16(z));
        end
    else %if sub is empty
        dsObj(i).subs = sub;
    end
    
end

save([MPN 'dsObj.mat'],'dsObj','-v7.3')
clear vastSubs










