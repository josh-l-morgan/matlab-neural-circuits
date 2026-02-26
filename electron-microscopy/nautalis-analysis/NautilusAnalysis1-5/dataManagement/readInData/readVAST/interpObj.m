function[dsObj] = interpObj(MPN,dSamp)

%% downSample object by 8

%MPN = 'D:\LGNs1\Segmentation\VAST\S8\joshm\export+14+04+27_mat\'

disp('loading vastSubs')
load([MPN 'vastSubs.mat'])
if ~exist('dSamp','var')
dSamp = [1 1 1];
end

%%

scale = 1./dSamp;
c = 0;
for y = 0:.5:ceil(scale(1))
    for x = .0:.5:ceil(scale(2));
        for z = 0:.5:ceil(scale(3));
            c = c+1;
           shiftScale(c,:) = [y x z];
        end
    end
end



%% down sample



for i = 1:length(vastSubs)
    
    disp(sprintf('%d of %d',i,length(vastSubs)));
    tic
    sub = double(vastSubs{i});
    if ~isempty(sub)
        sub(:,1) = round(double(vastSubs{i}(:,1))/dSamp(1));
        sub(:,2) = round(double(vastSubs{i}(:,2))/dSamp(2));
        sub(:,3) = round(double(vastSubs{i}(:,3))/dSamp(3));
        
        subSize = size(sub,1);
        newSub = zeros(size(sub,1)*size(shiftScale,1),3);
        
        
        
        
        sub(sub<1) = 1; %%!!! Problem getting subs should never get zeros
        maxSub = max(sub,[],1);
        if ~isempty(sub)
            inds = sub2ind(maxSub,sub(:,1),sub(:,2),sub(:,3));
            uInds = unique(inds);
            toc
            tic
            if length(uInds)>1
               hInds = hist(inds,uInds);
               
            else
                hInds = length(inds)
            end
            toc
            tic
            [y x z] = ind2sub(maxSub,uInds);
            dsObj(i).subs = cat(2,uint16(y),uint16(x),uint16(z));
            %dsObj(i).n = uint16(hInds);
            toc
        end
    else %if sub is empty
        dsObj(i).subs = sub;
    end
    
end

save([MPN 'dsObj.mat'],'dsObj','-v7.3')
clear vastSubs
