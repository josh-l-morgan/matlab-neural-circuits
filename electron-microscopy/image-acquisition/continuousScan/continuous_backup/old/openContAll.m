%% Change MAG


if ~exist('FibicsOn','var')
    TPN = GetMyDir;
    %% Activate ActiveX
    sm = actxserver('VBComObjectWrapperForZeissAPI.KHZeissSEMWrapperComClass')
    sm.InitialiseRemoting
    sm.Set_PassedTypeSingle('AP_MAG',25);
    sm.Fibics_Initialise();
    sprintf(' Fibics Initializing, pausing 15 seconds...')
    pause(15)
    FibicsOn = 1;
end
%%
startTime = clock;
Xstart =      109.9999/1000;
Xstop =     11.6591 / 1000; 

FOV = 4500;
W = 15000;
H = 15000;
dwell = .1;
WriteBase = [TPN 'test'];
Ystop =      4.4261/1000;%sm.Get_ReturnTypeSingle('AP_STAGE_GOTO_Y')
Ystart =      109.9999/1000;
Xstep = FOV / -1000000;
Xpos = Xstart: Xstep :Xstop;

for x = 1:length(Xpos)

WriteTo = [WriteBase num2str(x) '.tif'];

sm.Set_PassedTypeSingle('AP_STAGE_GOTO_X',Xpos(x));
smwait(sm,'DP_STAGE_IS');

sm.Set_PassedTypeSingle('AP_STAGE_GOTO_Y',Ystart);
smwait(sm,'DP_STAGE_IS');
'hi'
sm.Fibics_WriteFOV(FOV)
'start image'
sm.Fibics_AcquireImage(W,H,dwell,WriteTo);
'image going'
pause(.01)
sm.Set_PassedTypeSingle('AP_STAGE_GOTO_Y',Ystop);

   while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_STAGE_IS'));
    sm.Get_ReturnTypeSingle('AP_STAGE_AT_X')*1000;
pause(.01)
   end

while(sm.Fibics_IsBusy),  pause(.2),  end
%Yend= sm.Get_ReturnTypeSingle('AP_STAGE_AT_Y')*1000
sprintf('imaged strip %d of %d',x,length(Xpos))
end

    sm.Get_ReturnTypeSingle('AP_STAGE_AT_Y')*1000

stopTime= clock;

duration = stopTime - startTime

%% Build Strips

nams = getPics(TPN);
colormap gray(255)
%isoDir = [TPN(1:end-1) 'build'];
%mkdir(isoDir)
isoDir = [TPN(1:end - 1) 'iso\'];
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










