TPN = GetMyDir;
nams = getPics(TPN);
colormap gray(255)
%isoDir = [TPN(1:end-1) 'build'];
%mkdir(isoDir)
isoDir = GetMyDir;
xFilt = gaus3d([1 50 1],20);

for i = 1:length(nams)

    % I = imread([TPN nams{i}]);
    nam = nams{i};
    info = imfinfo([TPN nam]);
    width = info.Width;
    Cols = [1: 50 : width(1) width(1)];
    Rows = info.Height;
    I = 255-imread([TPN nams{i}],'PixelRegion',{[1,1,Rows(1)],[1,50,width(1)]});
    imwrite(I,[isoDir '\' nam])

    i
end



%%
buildDir = [isoDir(1:end -1) 'build'];
if ~exist(buildDir), mkdir(buildDir), end
nams = getPics(isoDir);
info = imfinfo([isoDir nams{1}]);

height = info.Height;
width = info.Width;
iNum = length(nams);
bI = zeros(height(1), width * iNum,'uint8');
Is = zeros(height(1), width, iNum,'uint8');
for i = 1:length(nams)
    nam = nams{i};
    I = imread([isoDir nam]);
    Is(:,:,i) = I;
end

%%
shiftRange = [0 : 1000];

for i = 2:size(Is,3)
    posProd = shiftRange * 0;
    negProd = shiftRange * 0;
    ref = double(Is(:,end,i-1));
    posSlide = double(Is(:,1,i));
    negSlide = posSlide; maxNeg = 0; maxProd = 0;
    for s = 1:1000
        posSlide = circshift(posSlide,[1 0]);
        negSlide = circshift(negSlide,[-1  0]);
        posProd(s) = sum(posSlide .* ref);
        maxProd = max(maxProd,posProd(s));

        negProd(s) = sum(negSlide .* ref);
        maxNeg = max(maxNeg,negProd(s));


        %image(posSlide),pause(.01)
        subplot(1,1,1)
        %         image([posSlide ref negSlide])
        %         subplot(1,2,2)
        plot(posProd,'r'),hold on
        plot(negProd,'g'),hold off
        pause(.01)
    end


    maxPos = max(posProd);
    maxNeg = max(negProd);
    if maxPos > maxNeg
        bestShift = find(posProd == maxPos);
    else
        bestShift = find(negProd == maxNeg) * -1;
    end
    bestShift

    %%
    Is(:,:,i) = circshift(Is(:,:,i),[bestShift 0]);
end

%% build
for i = 1:size(Is,3)
    bI( : , width * ( i - 1) + 1 : width * i) = Is(:,:,i);
end

image(bI)
if ~exist(buildDir),mkdir(buildDir),end
imwrite(bI,[buildDir '\bI.tif'])









