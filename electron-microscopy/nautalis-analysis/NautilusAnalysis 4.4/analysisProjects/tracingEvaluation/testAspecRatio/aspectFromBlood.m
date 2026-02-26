
clear all
sourceDir = 'D:\LGNs1\PSC_alignments\joshmProcess\S32intermediateSingleList_ds_medZ_filtBlood4_filtCBV3d_colorCBV\';

pngName = dir([sourceDir '*.tif']);

for p = 1:length(pngName)
    nam = pngName(p).name;
    colI(:,:,:,p) = imread([sourceDir nam]);
    
end

bvI = squeeze(colI(:,:,1,:));

image(sum(bvI,3)*.03)
image(squeeze(sum(bvI,2))*.03)

%%
binWidth = 50;
c = 0;
medC = [];
for y = 1:binWidth:size(bvI,1)-binWidth;
    y
    for x = 1:binWidth:size(bvI,2)-binWidth
        for z = 1:binWidth:size(bvI,3)-binWidth
            
            disp(sprintf('%d %d %d',y,x,z))
            grab = bvI(y:y+binWidth-1,x:x+binWidth-1,z:z+binWidth-1);
            
            [yg xg zg] = ind2sub(size(grab),find(grab>0));
            if ~(isempty(yg))
                
                
                
                [labGrab num] = bwlabeln(grab);
                for n = 1:num
                    c = c+1;
                    [yg xg zg] = ind2sub(size(labGrab),find(labGrab>0));
                    medC(c,:) = [median(yg)+y-1 median(xg)+x-1 median(zg)+z-1];
                    
                end
                
               
            end
        end
    end
end



%% 

testY = [3:.1:8];
testX = [4];
testZ = [20:40];

[by bx bz] = ind2sub(size(bvI),find(bvI>0));
maxDist = 20;
clear latents
latMat = zeros(length(testY), length(testZ), size(medC,1),3);
for i = 1:size(medC,1)
    disp(sprintf('%d of %d',i,size(medC,1)));
    dists = sqrt((by-medC(i,1)).^2 + (bx-medC(i,2)).^2 + ...
        (bz-medC(i,3)).^2 );
    useC = [by(dists<maxDist,:)  bx(dists<maxDist,:)  ...
        bz(dists<maxDist,:)];
    
    
    
    t = 0;
    for ty = 1:length(testY)
        for tx = 1:length(testX)
            for tz = 1:length(testZ)
                    testC = [useC(:,1)*testY(ty) useC(:,2)*testX(tx) useC(:,3)*testZ(tz)];
                    [coef, score, latent] = princomp(testC);
                    t = t+1;
              %latents(i,t,:) = latent;     
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





