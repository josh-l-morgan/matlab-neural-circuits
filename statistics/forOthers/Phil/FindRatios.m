

%%Bad image
SPN = 'Y:\Active\sean\TwitchData\Images\KCNG\Processed\904-0 Processed\export1\'
%%OK image
SPN = 'Y:\Active\sean\TwitchData\Images\KCNG\Processed\904-0 Processed\export2\'
%%OK image
SPN = 'Y:\Active\sean\TwitchData\Images\KCNG\Processed\904-0 Processed\export3\'


shouldFilt = 1;


dSPN = dir([SPN '*.tif']);
nams = {dSPN.name};

clear I1 I2
for i = 1:length(nams)
    nam = nams{i};
    cs = regexp(nam,'_c');
    dot = regexp(nam,'.tif');
    col = str2num(nam(cs(1)+2:dot(1)-1));
    zs = regexp(nam,'_z');
    z = str2num(nam(zs(1)+2:cs(1)-1));
        
    I = double(imread([ SPN nams{i}]));
    
    if shouldFilt
        %         I = medfilt2(I,[3 3],'symmetric');
        %se = fspecial('gaussian',10,5
        %I = imfilt
        I = imgaussfilt(I,10);
    end
    
    
    if col == 1
        I1(:,:,z) = I;
    elseif col == 2
        I2(:,:,z) = I;
    end
        
    
end
  
I1 = I1 - min(I1(:));
I1 = I1/mean(I1(:));
I2 = I2 - min(I2(:));
I2 = I2/mean(I2(:));

subplot(1,1,1)
for i = 1:size(I1,3)
    Ic = cat(3,I1(:,:,i),I2(:,:,i),I2(:,:,i));
    image(uint8(Ic*20))
    pause(.1)
end

%%%%%%%%%%%%%%%%%
val1 = I1(:);
val1 = val1-min(val1(:));
val1 = val1/mean(val1);
val2 = I2(:);
val2 = val2-min(val2(:));
val2 = val2/mean(val2);

valM = val1+val2;
[valM idx] = sort(valM,'ascend');
val2 = val2(idx);
val1 = val1(idx);

rat = val1./val2;

subplot(3,1,1)
plot(val1,'r')
subplot(3,1,2)
plot(val2,'g')
subplot(3,1,3)
plot(rat,'b')

allMin = min([val1 val2]);
allMax = max([val1 val2]);
hRange = [allMin :.1: allMax];

hVal1 = hist(val1,hRange);
hVal2 = hist(val2,hRange);

subplot(1,1,1)
plot(hRange,hVal1,'r');
hold on
plot(hRange,hVal2),'g';
hold off

%% Divide
binNum = 10;
numV = length(val1);
vBins = ceil([1 [1:binNum]* numV/binNum]);
clear meanRat meanVal1 meanVal2 varRat
for i = 1:length(vBins)-1
    meanRat(i) = mean(rat(vBins(i):vBins(i+1)));  
    varRat(i) = var(rat(vBins(i):vBins(i+1)));  
    meanVal1(i) = mean(val1(vBins(i):vBins(i+1)));  
    meanVal2(i) = mean(val2(vBins(i):vBins(i+1)));  
end

subplot(3,1,1)
plot(meanVal1,'r')
hold on
plot(meanVal2,'g')
subplot(3,1,2)
plot(meanRat)
hold on
plot(meanVal1./meanVal2,'r')
hold off
subplot(3,1,3)
plot(varRat./meanRat)
%ylim([0 3])

%% Optimize
clf
b1 = [-1:.2:3];
b2 = [-1:.2:1];
c = 0;
useBin = [5:10];
clear eVar iVar
for y = 1:length(b1)
    for x = 1:length(b2)
        
        clear nMeanRat nMeanVal1 nMeanVal2 nVarRat
        
        nVal1 = val1 + b1(y);
        
        nVal2 = val2 + b2(x);
        nRat = nVal1./nVal2;
        for i = 1:length(vBins)-1
            nMeanRat(i) = mean(nRat(vBins(i):vBins(i+1)));
            nVarRat(i) = var(nRat(vBins(i):vBins(i+1)));
            nMeanVal1(i) = mean(nVal1(vBins(i):vBins(i+1)));
            nMeanVal2(i) = mean(nVal2(vBins(i):vBins(i+1)));
        end
        
        iVar(y,x) = mean(nVarRat(useBin));
        eVar(y,x) = var(nMeanRat(useBin));
        
        
%         c = c+1;
%         subplot(1,1,1)
%         %subplot(11,11,c)
%         plot(nMeanRat)
%         hold off
%         %plot(nVarRat./nMeanRat,'r')
%         ylim([-10 10])
%         pause(.01)
        
        
    end
end
hold off

subplot(2,1,1)
colormap jet(100)
image(iVar/max(iVar(:))*100)
subplot(2,1,2)
image(eVar/max(eVar(:))*100)




