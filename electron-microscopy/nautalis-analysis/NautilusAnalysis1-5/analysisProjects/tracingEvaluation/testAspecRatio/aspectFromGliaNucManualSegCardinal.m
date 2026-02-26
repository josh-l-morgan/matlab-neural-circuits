
clear all
load('MPN.mat')
load([MPN 'obI.mat'])
load([MPN 'vastSubs.mat']);
names = obI.colStruc.names;
bvOb = [];
for i = 1:length(names)
    nam = names{i};
    if sum(regexp(nam,'glnuc'))
        bvOb = [bvOb i];
    end
end

clear lengthOb
for i = 1:length(bvOb)
    useC = vastSubs{bvOb(i)};
    lengthOb(i) = size(useC,1);
end
bvOb = bvOb(lengthOb>0);


%%

testY = [3:.1:8];
testX = [4];
testZ = [30];

latMat = zeros(length(testY), length(testZ),length(bvOb),3);
for i = 1: length(bvOb)
    disp(sprintf('%d of %d',i,length(bvOb)));
    
    
    useC = vastSubs{bvOb(i)};
    for d = 1:3
        useC(:,d) = useC(:,d)- min(useC(:,d))+1;
    end
    maxC = max(useC,[],1);
    
    if 0
        useV = zeros(maxC);
        inds = sub2ind(maxC,useC(:,1), useC(:,2),useC(:,3));
        useV(inds) = 1;
        for p = 1:maxC(3);
            tempI = useV(:,:,p);
            SE = strel('disk',3);
            tempI2 = imopen(tempI,SE);
            useV(:,:,p) = tempI2;
            
            
            subplot(2,1,1)
            image(tempI*100)
            subplot(2,1,2)
            image(tempI2*100)
            pause(.01)
            
        end
        newInds = find(useV>0);
        [ys xs zs] = ind2sub(maxC,newInds);
        newC = [ys xs zs];
    else
        newC = useC;
    end
    
    scaleC = scaleSubs(newC,[1 1 1/4]);
    downSamp = 4;
    smallSub = shrinkSub(scaleC,downSamp);
    
    fv = subVolFV(scaleC,[],1);
    %      pause
    
    tic
    t = 0;
    for tx = 1:length(testX)
        for tz = 1:length(testZ)
            for ty = 1:length(testY)
                
                testC = [useC(:,1)*testY(ty) useC(:,2)*testX(tx) useC(:,3)*testZ(tz)];
                downSamp = 2;
                smallSub = shrinkSub(testC,downSamp);
                %[coef, score, latent] = princomp(smallSub);
                t = t+1;
                %latents(i,t,:) = latent;
                latMat(ty,tz,i,:) = var(smalSub);
            end
        end
    end
    toc
    
    
end


isZer = [];
for i = 1:size(latMat,3);
    samp = latMat(:,:,i,:);
    isZer(i) = sum(sum(sum(samp==0)));
end
latMat = latMat(:,:,isZer ==0,:);

latRat = (latMat(:,:,:,2)-latMat(:,:,:,3))./latMat(:,:,:,2);

rat1 = abs((latMat(:,:,:,1)-latMat(:,:,:,2))./latMat(:,:,:,2))
rat2 = abs((latMat(:,:,:,3)-latMat(:,:,:,2))./latMat(:,:,:,2))


meanLatRat = mean(latRat,3);
showRat = meanLatRat - min(meanLatRat(:));

image(showRat * 100/max(showRat(:)))

[miny minz] = find(meanLatRat == min(meanLatRat(:)));

bestY = testY(miny)
bestZ = testZ(minz)





