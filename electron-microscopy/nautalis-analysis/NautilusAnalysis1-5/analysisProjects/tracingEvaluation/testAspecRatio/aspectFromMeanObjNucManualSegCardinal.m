
clear all
load('MPN.mat')
load([MPN 'obI.mat'])
load([MPN 'dsObj.mat']);
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
    useC = dsObj(bvOb(i)).subs;
    lengthOb(i) = size(useC,1);
end
bvOb = bvOb(lengthOb>0);


%%
defaultAspect = [ 1 1 1];

testY = [.5:.1:1.5];
testX = [1];
testZ = [.5:.1:1.5];


latMat = zeros(length(testY), length(testZ),3);
catSub = [];
for i = 1: length(bvOb)
    disp(sprintf('%d of %d',i,length(bvOb)));
    
    
    useC = dsObj(bvOb(i)).subs;
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
    
    
    %fv = subVolFV(round(newC-min(newC(:))+1),[],1);
    %      pause
    
    for ty = 1:length(testY)
        testC = [newC(:,1)*testY(ty) newC(:,2)* testX(1) newC(:,3) * defaultAspect(3)];
        downSamp = 2;
        smallSub = shrinkSub(testC,downSamp);
        normSub = [smallSub(:,1) - mean(smallSub(:,1)) ...
            smallSub(:,2) - mean(smallSub(:,2)) smallSub(:,3) - mean(smallSub(:,3)) ];
        
        try
            catSubY{ty} =  cat(1,catSubY{ty} , normSub);
        catch err
            catSubY{ty}  = normSub;
        end
    end
    
    for tz = 1:length(testZ)
        
        testC = [newC(:,1)*defaultAspect(1) newC(:,2)* testX(1) newC(:,3)*testZ(tz)];
        downSamp = 2;
        smallSub = shrinkSub(testC,downSamp);
        normSub = [smallSub(:,1) - mean(smallSub(:,1)) ...
            smallSub(:,2) - mean(smallSub(:,2)) smallSub(:,3) - mean(smallSub(:,3)) ];
        
        try   catSubZ{tz} =  cat(1,catSub{tz} , normSub);
        catch err
            catSubZ{tz}  = normSub;
        end
    end
    toc
    
end

latMatY = zeros(length(testY),3);
for ty = 1:length(testY)
    
    testSub = catSubY{ty};
    latMatY(ty,:) = var(double(testSub));
end

latMatZ = zeros(length(testZ),3);
for tz = 1:length(testZ)
    
    testSub = catSubZ{tz};
    latMatZ(tz,:) = var(double(testSub));
end


%
% isZer = [];
% for i = 1:size(latMat,3);
%     samp = latMat(:,:,i,:);
%     isZer(i) = sum(sum(sum(samp==0)));
% end
% latMat = latMat(:,:,isZer ==0,:);


rat1 = abs((latMatY(:,1)-latMatY(:,2))./latMatY(:,2))
rat2 = abs((latMatZ(:,3)-latMatZ(:,2))./latMatZ(:,2))



bestY = testY(find(rat1 == min(rat1)))
bestZ = testZ(find(rat2 == min(rat2)))





