
clear all
load('MPN.mat')
load([MPN 'obI.mat'])
load([MPN 'dsObj.mat']);
names = obI.colStruc.names;
bvOb = [];
for i = 1:length(names)
    nam = names{i};
    if sum(regexp(nam,'bldv'))
        bvOb = [bvOb i];
    end
end

for i = 1:length(bvOb)
        useC = dsObj(bvOb(i)).subs;
        lengthOb(i) = size(useC,1);
end
bvOb = bvOb(lengthOb>0);


%% 

testY = [.5:.1:1.5];
testX = [1];
testZ = [1];

latMat = zeros(length(testY), length(testZ),length(bvOb),3);
for i = 1:length(bvOb)
    disp(sprintf('%d of %d',i,length(bvOb)));
    
    
    useC = double(dsObj(bvOb(i)).subs);
    for d = 1:3
        useC(:,d) = useC(:,d)- min(useC(:,d))+1;
    end
    maxC = max(useC,[],1);
    useV = zeros(maxC);
    inds = sub2ind(maxC,useC(:,1), useC(:,2),useC(:,3));
    useV(inds) = 1;
    for p = 1:size(useV,3);
       tempI = useV(:,:,p);
       SE = strel('disk',3);
       tempI2 = imopen(tempI,SE);
       useV(:,:,p) = tempI2; 
       
%        
%        subplot(2,1,1)
%        image(tempI*100)
%        subplot(2,1,2)
%        image(tempI2*100)
%        pause(.01)
       
    end
    newInds = find(useV>0);
    [ys xs zs] = ind2sub(maxC,newInds);
    newC = [ys xs zs];
     image(squeeze(sum(useV,2))*3)
    
    downSamp = 4;
    smallSub = shrinkSub(newC,downSamp);

     fv = subVolFV(smallSub,[],1);
     pause(.01)
    
    t = 0;
    for ty = 1:length(testY)
        for tx = 1:length(testX)
            parfor tz = 1:length(testZ)
                    testC = [newC(:,1)*testY(ty) newC(:,2)*testX(tx) newC(:,3)*testZ(tz)];
                     downSamp = 2;
                    smallSub = shrinkSub(testC,downSamp);
                    [coef, score, latent] = princomp(smallSub);
              latMat(ty,tz,i,:) = latent;
            end
        end
    end
    
    
    
end


isZer = [];
for i = 1:size(latMat,3);
    
    samp = latMat(:,:,i,:);
    isZer(i) = sum(sum(sum(samp==0)));
    
    
end
latMat = latMat(:,:,isZer ==0,:);

latRat = (latMat(:,:,:,2)-latMat(:,:,:,3))./latMat(:,:,:,2);

meanLatRat = mean(latRat,3);
showRat = meanLatRat - min(meanLatRat(:));

image(showRat * 100/max(showRat(:)))

[miny minz] = find(meanLatRat == min(meanLatRat(:)));

bestY = testY(miny)
bestZ = testZ(minz)





