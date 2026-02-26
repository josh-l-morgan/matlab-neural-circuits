colormap gray(256)

[TFN TPN] = GetMyFile;


Iraw = imread([TPN TFN]);
I = Iraw;
I = I(500:1000,400:900);
I(1:10,:) = 0;
I(end -10:end,:) = 0;
I(:,1:10) = 0;
I(:,end-10:end) = 0;

subplot(1,2,1)

image(I)
subplot(1,2,2)
%I = 255-I;
[ys, xs] = size(I);

%% Filter for membranes
gKern1 = gaus3d([15 15 1],1);
c1 = fastCon(I,gKern1);
image(fitH(c1))

gKern2 = gaus3d([200 200 1],30);
c2 = fastCon(I,gKern2);
image(fitH(c2))

difI = c1 - c2;
image(fitH(difI))
%% Edge
gKern1 = gaus3d([15 15 1],1);
c1 = fastCon(I,gKern1);
image(fitH(c1))
eI = edge(c1,'Canny',[.000005 .3],1);
image(fitH(eI*1000+c1))

%% close dark
% eI(1:10,:) = 0;
% eI(end -10:end,:) = 0;
% eI(:,1:10) = 0;
% eI(:,end-10:end) = 0;
cI = eI;
sRad = 10;
for r = 1:1
    pI = bwperim(cI);
    edgeIDs = find(cI);
    diskE = gaus3d([sRad * 2 + 1,sRad * 2 + 1,1],3);
    diskE = diskE>=diskE(sRad,1);
    [y x] = find(diskE);
    yx = [y x]-sRad;
    for i = 1:length(edgeIDs)
        dVal = dI(edgeIDs(i));
        [cy cx] = ind2sub(size(I),edgeIDs(i));
        y = yx(:,1) + cy; x = yx(:,2) + cx;
        uyx = wall([y x],size(I));
        iyx = sub2ind(size(I),uyx(:,1),uyx(:,2));
        sVal = dI(iyx);  % get values of surround
        cI(iyx(sVal>dVal)) = 1;   %set dark surround to 1
        %image(fitH(cI + eI)),pause
    end
    image(fitH(cI + eI)),pause(1)
end
%%
gKern2 = gaus3d([15 15 1],5);
c2 = fastCon(I,gKern1);
image(fitH(c2))
c2 = max(c2(:)) - c2;
image(fitH(c2))
%dI2 = imhmin(dI,20);

wI = watershed(dI2,8);
image((wI==0)*1000 + dI2 * .8);
pI = I * 0;
pI((wI==0) & (cI>0)) = 1;
image(pI * 100)
%Arbitrary trimming
% image((wI==0) * 100 + (xI)*100)
% xI(wI==0) = 0;

%% Label

lI = bwlabel(~pI,4);
colI = uint8(cat(3,I,I,I))/2;
myCol = hsv(100)*100;
skipCol = 5 + round(rand*10);
for i = 1:max(lI(:))
    [y x] = find(lI==i);
    uCol = myCol(mod(i,size(myCol))+1,:);
    uCol = myCol(mod(i*skipCol,size(myCol,1))+1,:);
    for c = 1:3
    colI(:,:,c) = colI(:,:,c)+ uint8(lI==i)*uCol(c);
    end

end
image(colI),pause(.1)







