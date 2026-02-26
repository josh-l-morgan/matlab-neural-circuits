
clear all

for condition = 1:3
    clear rec
    
    if  condition == 1 % P32 LGNs1
        
        SPN = 'D:\LGN\2DBouton_Segmentations\2DboutonCharacterizations\P7_GxK\'
        voxSize = [4 4 40];        
        backgroundID = 0;
        rgcID = 1;
        areaID = 4;
        cbID = 16;
        unkID = 2;
        bvID = 7;
                
    elseif condition == 2
        
        SPN = 'D:\LGN\2DBouton_Segmentations\2DboutonCharacterizations\P11_Ixd_W1_segImage\'
        voxSize = [4 4 40];        
        backgroundID = 0;
        rgcID = 1;
        areaID = 4;
        cbID = 16;
        unkID = 2;
        bvID = 7;
                
    elseif condition == 3
        
        SPN = 'D:\LGN\2DBouton_Segmentations\2DboutonCharacterizations\P32_LGNs1_segImage\'
        voxSize = [16 16 30];        
        backgroundID = 0;
        rgcID = 1;
        areaID = 4;
        cbID = 16;
        unkID = 2;
        bvID = 5;
                
    end
    
    
    searchRad = 3; %radius in um
    
    
    Idir = dir([SPN '*.png']);
    Inam = {Idir.name};
    
    pixArea = voxSize(1)/1000 * voxSize(2) / 1000;
    
    pixRad = round(searchRad / voxSize(1) * 1000);
    
    disk = ones(pixRad*2+1,pixRad*2+1);
    [y x] = find(disk);
    dists = sqrt((y-mean(y)).^2 + (x-mean(x)).^2);
    disk(dists>pixRad) = 0;
    diskY = y(dists<=pixRad)-pixRad-1;
    diskX =  x(dists<=pixRad)-pixRad-1;
    
    for i = 1:length(Inam)
        nam = Inam{i};
        I = imread([SPN nam]);
        
        totAll = sum(I(:)>0);
        [y x] = find(I>0);
        centAll = [mean(y) mean(x)];
        
        %%get Area edges
        areaI = I>0;
        [y x] = find(areaI);
        cent = [mean(y) mean(x)];
        
        rgcI = I == rgcID;
        totRGC = sum(rgcI(:));
        totNeuroPil = sum(I(:) == areaID);
        totCB = sum(I(:) == cbID);
        totBV = sum(I(:) == bvID);
        
        
        props = regionprops(rgcI,'Area','centroid');
        cents = cat(1,props.Centroid);
        numRGC = length(props);
        
        
        percentRGC = totRGC/totAll * 100;
        percentRGCneuroPil = totRGC/totNeuroPil * 100;
        areaAll = totAll * pixArea;
        areaRGC = totRGC * pixArea;
        areaNeuroPil = totNeuroPil * pixArea;
        areaCB = totCB * pixArea;
        areaBV = totCB * pixArea;
        
        RGCper100UM = numRGC / areaAll * 100;
        RGCper100UMneuroPil = numRGC / areaNeuroPil * 100;
        
        res.numRGC = numRGC;
        res.areaAll = areaAll;
        res.areaNeuropil = areaNeuroPil;
        res.percentNeuropil = areaNeuroPil/areaAll;
        res.RGCper100um = RGCper100UM;
        res.percentRGC = percentRGC;
        res.RGCper100umNeuropil = RGCper100UMneuroPil;
        res.percentRGCneuroPil = percentRGCneuroPil;
        res.rgcAreas = cat(1,props.Area)* pixArea;
        res.rgcCents = cat(1,props.Centroid);
        diskI = I*20;
        
        centUM = [cents(:,1) * voxSize(1)/1000 ...
            cents(:,2) * voxSize(2)/1000];
        clear bout
        for b = 1: numRGC
            
            shiftDisk = round(sub2ind(size(I),...
                diskY + round(cents(b,2)),diskX + round(cents(b,1))));
            diskI(shiftDisk) = rand*100;
            
            valI = I(shiftDisk);
            
            bout(b).cusioned = (sum(valI==0))/length(valI)< 0.05;
            bout(b).percentRGC = sum(valI==rgcID)/sum(valI>0)*100; %percent of nearby region that is RGC 
            dists = sqrt((centUM(:,1)-centUM(b,1)).^2 + ...
                (centUM(:,2)-centUM(b,2)).^2);
            bout(b).near = sum(dists<=searchRad)-1;
        end
        
        
        image(diskI),pause(.01)
        pause(.01)
        cusionedBouts = ([bout.cusioned])>0;
        nearPercentRGCs = [bout(cusionedBouts).percentRGC];
        nearNums = [bout(cusionedBouts).near];
        nearDens = nearNums/ (pi * searchRad^2) * 100;
        
        nearHistRange = [0:20];
        nearHist = hist(nearNums,nearHistRange);
        bar(nearHist)
        
        quant.meanRGCarea = mean(res.rgcAreas);
        quant.medianRGCarea = median(res.rgcAreas);
        quant.varRGCarea = var(res.rgcAreas);
        quant.seRGCarea = var(res.rgcAreas)/sqrt(length(res.rgcAreas));
        
        
        
        boutCus = cat(1,bout.cusioned);
        boutNear = cat(1,bout.near);
        boutPercentRGC = cat(1,bout.percentRGC);
        
        quant.boutNear = boutNear(boutCus);
        quant.percentRGC = boutPercentRGC(boutCus);
        
        quant.meanRGCnear = mean(quant.boutNear);
        quant.medianRGCnear  = median(quant.boutNear);
        quant.meanNearPercentRGC = mean(quant.percentRGC);
        quant.medianNearPercentRGC = median(quant.percentRGC);
        
        
        
        scatter(quant.boutNear,quant.percentRGC)
            
    end
    
    r(condition).res = res;
    r(condition).bout = bout;
    r(condition).quant = quant;
    pause(.01)
end


%%
  
areaRanksum = ranksum(r(2).res.rgcAreas, r(3).res.rgcAreas)

vals = r(2).res.rgcAreas
val1 = abs(vals - mean(vals))/mean(vals)
vals = r(3).res.rgcAreas
val2 = abs(vals - mean(vals))/mean(vals)
ranksum(val1,val2)


vals = r(2).quant.boutNear
mean(vals)
length(vals)
std(vals/sqrt(length(vals)))



nearsRanksum = ranksum(r(2).quant.boutNear, r(3).quant.boutNear) 

areaRange = 0:.3:4;
histArea1 = hist(r(2).res.rgcAreas,areaRange);
histArea1 = histArea1/sum(histArea1);
histArea2 = hist(r(3).res.rgcAreas,areaRange);
histArea2 = histArea2/sum(histArea2);

bar(areaRange, [histArea2' histArea1'],1.3)



nearRange = 0:2:20;
histNear1 = hist(r(2).quant.boutNear,nearRange);
histNear1 = histNear1/sum(histNear1);
histNear2 = hist(r(3).quant.boutNear,nearRange);
histNear2 = histNear2/sum(histNear2);

bar(nearRange, [histNear2' histNear1'],1.3)









