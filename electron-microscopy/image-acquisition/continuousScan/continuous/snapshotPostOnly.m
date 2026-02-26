%% Build map out of acquired images

TPN = GetMyDir
rawDir = [TPN 'raw\'];

nams = getPics(rawDir);
colormap gray(255)
%isoDir = [TPN(1:end-1) 'build'];
%mkdir(isoDir)
isoDir = [TPN 'isoDir\'];
if ~exist(isoDir),mkdir(isoDir),end


for i = 1:length(nams)
    sprintf('reading %d of %d',i,length(nams))

    nam = nams{i};
    info = imfinfo([rawDir nam]);
    width = info.Width;
    Cols = [1: 31 : width(1) width(1)];
    Rows = info.Height;
    Is{i} = 255-imread([rawDir nams{i}],'PixelRegion',{[1,1,Rows(1)],[1,31,width(1)]});

end

%% Correct flips

maxdim = [0 0];

%shrinkFlip = Ymove/(Ymove + FOV * .75);
shrinkFlip = 1.0712;
for i = 1: length(Is)
    flip = 1 - mod(i,2);
    I = Is{i};
    if flip
        I = flipud(I);
        I = circshift(I,[-2600 0]); 
        I = imresize(I,[size(I,1)*shrinkFlip,size(I,2)],'nearest');
      
    else
    end
    maxdim = max(maxdim,[size(I,1) size(I,2)]);
    %imwrite(I,[isoDir nam])
    image(I),pause(.01)
    It{i} = I;
end

%%
for i = 1 :length(It)
    subplot(1,length(It),i)
    image(It{i})
    YLim([0 maxdim(1)])
end

pause(1)
%% Align Strips
shiftRange = [0 : 3000];
Islide = It{1};

for i = 2:length(Is)
    posProd = shiftRange * 0;
    negProd = shiftRange * 0;
    Iref = Islide;
    ref = double(Iref(:,end));
    ref(maxdim(1)) = 0;
    Islide = It{i};
    posSlide = double(Islide(:,1));
    posSlide(maxdim(1)) = 0;
    negSlide = posSlide; 
    
    maxNeg = 0; maxProd = 0;
    sprintf('sliding %d of %d', i,size(Is,3))
    goN = 1; goP =1;
    for s = shiftRange
        if goP
            posSlide = circshift(posSlide,[1 0]);
            posProd(s+1) = sum(posSlide .* ref);
            maxProd = max(maxProd,posProd(s+1));
        end
        if goN
            negSlide = circshift(negSlide,[-1  0]);
            negProd(s+1) = sum(negSlide .* ref);
            maxNeg = max(maxNeg,negProd(s+1));
        end
        
        if s>200
            if mean((negProd(s-20:s-1)-negProd(s-19:s))>0)> 0.9
                goN = 0;
            end
            if mean((posProd(s-20:s-1)-posProd(s-19:s))>0)> 0.9
                goP = 0;
            end
            if ~goN & ~goP
                break
            end
        end
%         %image(posSlide),pause(.01)
%         subplot(1,2,1)
%         image([posSlide ref negSlide])
%         subplot(1,2,2)
%         plot(posProd,'r'),hold on
%         plot(negProd,'g'),hold off
%         pause(.01)
    end
    
    
    maxPos = max(posProd);
    maxNeg = max(negProd);
    if maxPos > maxNeg
        bestShift = find(posProd == maxPos)-1;
    else
        bestShift = find(negProd == maxNeg) * -1 - 1;
    end
    bestShift
    
    Is{i} = circshift(Islide,[bestShift 0]); % excecute shift
end

%% Build single image
subplot(1,1,1)
buildDir = [TPN 'build'];
if ~exist(buildDir), mkdir(buildDir), end

bI = zeros(maxdim(1), maxdim(2)*length(Is),'uint8');

startStrip = 1;
for i = 1:length(It)
    I = It{i};
    stopStrip = startStrip + size(I,2)-1;
    bI( 1:size(I,1) , startStrip : stopStrip) = I;
    startStrip = stopStrip + 1
end

image(bI)
if ~exist(buildDir),mkdir(buildDir),end
imwrite(bI,[buildDir '\bI.tif'])









