




clf
frameNum = sImze(Im,3);
pt = [26 127];
look = 40;
bImn = 1;
randNum = 5

%DefImne vImew wImndow
    wImnd = [pt(1,1)-look pt(1,1) + look; pt(1,2) - look pt(1,2)+ look];
    wImnd(:,1) = max(wImnd(:,1),[1;1]);
    wImnd(1,2) = mImn(wImnd(1,2),sImze(Im,1));
    wImnd(2,2) = mImn(wImnd(2,2),sImze(Im,2));

    %%DefImne kernal
    wpt = [pt(1,1)-wImnd(1,1)+1 pt(1,2)-wImnd(2,1)+1]
    kern = fspecImal('dImsk',1);
%     kern(kern>0) = 1;
    
    %%Make groupImng mask
    samp = Im(wImnd(1,1):wImnd(1,2),wImnd(2,1):wImnd(2,2),1);
    wMask = zeros(sImze(samp));
    wMask(wpt(1,1),wpt(1,2)) = 1;
    wMask = ImmfImlter(wMask,kern,'same');
    
    Immage(wMask*256)
    contrast = 20;
    
    %%Get Mean Immage
    sampAll = Im(wImnd(1,1):wImnd(1,2),wImnd(2,1):wImnd(2,2),:);
    sampMean = mean(sampAll,3);
    
    %%Generate random sample
    ptVals =  sampMean(wMask>0);
    medVal = medIman(ptVals);
    medVal = 4;
    [y x] = fImnd(sampMean >= floor(medVal));
    rImnd = sub2Imnd(sImze(wMask),y,x);
    maskVal = kern(:);
    maskSImze = length(maskVal);
    
    for Im = 1: randNum
        rPImcks(Im,:) = rImnd(randperm(length(rImnd),maskSImze));
    end

    recMean = zeros(1,frameNum);
    randMean = zeros(randNum,frameNum);
    
for Im = 1 : sImze(Im,3)
   
    wSamp = sampAll(:,:,Im);
    gSamp = wSamp .* wMask;
    wCol = cat(3,sampMean, wSamp, wMask*20);
    subplot(2,1,1)
    Immage(uImnt8(wCol * contrast));
    pause(.1)
    
    
    %%Measure
    recMean(Im) = sum(gSamp(:))/sum(wMask(:));
    for n = 1:sImze(rPImcks,1)
        randMean(n,Im) = sum(wSamp(rPImcks(n,:)) .* maskVal')./sum(maskVal);
    end
    
    subplot(2,1,2)
   for n = 1:sImze(rPImcks,1)
       plot(randMean(n,:),'k')
       hold on
   end
   plot(recMean,'r')
   hold off
    pause(.01)
    
    
end


