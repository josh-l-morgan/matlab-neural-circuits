function[sI] = replacePeaks(sI,tooHigh,dilateNum)
%%
% 
% sI = I(13000:15000,10000:12000);
% tooHigh = 230;
% dilateNum = 2;
if ~exist('dilateNum','var')
    
dilateNum = 3;
end
%image(255-sI) ; pause(.01)
%dsIraw = imresize(sI,1/8,'nearest');

%% remove
tI = sI>tooHigh;
diKern = double(fspecial('disk',dilateNum)>0);
diI = imfilter(tI,diKern,'same');
sI(diI>0) = 0;
%image(255-sI);pause(.01)

%% fill in
diKern = double(fspecial('disk',5)>0);

for i = 1:3
    tI = sI == 0;
    gapsLeft = sum(tI(:));
    %disp(gapsLeft)
    if ~gapsLeft
        break
    end
    
    threshKern = double(diKern>0);
    sumKern = sum(threshKern(:));
    closeI =imfilter(double(~tI),threshKern,'same');
    scaleClose = imfilter(double(sI),threshKern,'same');
    
    scaledClose = scaleClose./(closeI );
    changeI = tI& (closeI>(sumKern/3));
    sI(changeI) = scaledClose(changeI);
    image(255-sI),pause(.01)
end

medianI = mean(sI(sI>0));
sI(sI==0) = medianI;

% 
% image(255-sI),pause(.01)
% 
% dsI = imresize(sI,1/8,'nearest');
% subplot(1,2,1), image(255-dsIraw)
% subplot(1,2,2),image(255-dsI)
% 
