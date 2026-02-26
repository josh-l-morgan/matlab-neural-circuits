%%testOverlaps


%%-Josh Morgan
[TPN] = GetMyDir;  %Get Tif
inams = getTifs(TPN);

for i = 1:inams - 1
    I1 = imread([TPN inams{i}]);
    lI1 = imread([TPN inams{i}]);
    I2 = imread([TPN inams{i+1}]);
    lI2 = imread([TPN inams{i+1}]);

    I1 = cat(3,I1,I1,I1); %Color image
    I2 = cat(3,I2,I2,I2);
    
    %%Make colormap
    myCol = hsv(max(lI(:))+1)*40;
    [r rix] = sort(rand(size(myCol,1)),1);
    myCol= myCol(rix,:);
    myCol = cat(1,[0 0 0],myCol);
    red = myCol(:,1); green = myCol(:,2); blue = myCol(:,3);
    skipCol = 5 + round(rand*10);

    col1 = uint8(cat(3,red(lI1+1),green(lI1+1),blue(lI1+1)));
    col2 = uint8(cat(3,red(lI2+1),green(lI2+1),blue(lI2+1)));
    
    while 1
    image(I1 * .7 + col1),pause(.01)
    
    pause
        image(I2 * .7 + col2),pause(.01)
    
    pause
    
    end
    
    
    [ys, xs] = size(I);
    subplot(1,2,1)
    image(I)
    subplot(1,2,2)

end