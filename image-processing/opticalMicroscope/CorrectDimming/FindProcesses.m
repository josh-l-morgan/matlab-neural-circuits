

SPN = 'H:\LxA\mx\TifStack\'
%SPN = 'H:\LxA\mx\redGreenGray\'
TPN = 'H:\LxA\mx\redGreenGray_tweak\'
if ~exist(TPN,'dir'),mkdir(TPN);end

dSPN = dir([SPN '*.tif']);
nams = {dSPN.name}

clear I
for i = length(nams):-1:1
    nam = nams{i};
    I(:,:,:,i) = double(imread([SPN nam]));
end


clear iMed
x = 1:size(I,4);
for i = 1:size(I,4)

    for c = 1:3
    Is = double(I(:,:,c,i));
    iMed(i,c) = median(Is(:));
    end

end
plot(iMed)

clf, hold on
pCol = [1 0 0; 0 1 0; 0 0 1];
for c = 1:3
    P = polyfit(x,iMed(:,c),1);
    yfit = P(1)*x+P(2);
    plot(x,iMed(:,c),'color',pCol(c,:))
    plot(x,yfit,'color',pCol(c,:))
    cScale(:,c) = max(yfit(:))./yfit;

end

I = I;
for i = 1:size(I,4)
    for c = 1:3
        I(:,:,c,i) = I(:,:,c,i) * cScale(i,c);
    end
end


clf
for i = 1:size(I,4)

    Is = I(:,:,:,i);
    image(uint8(Is))
    drawnow


end


%% Band pass filter

s1 = 2;
s2 = 10;
k1 = fspecial('gaussian',s1*4,s1);
k2 = fspecial('gaussian',s2*4,s2);


clear iMed
x = 1:size(I,4);
Ic1 = zeros(size(I,1),size(I,2), size(I,4),'single');
Ic2 = Ic1;
Ic3 = Ic1;

for i = 1:size(Ic1,3)
    disp(sprintf('bandpass filtering section %d of %d',i,size(Ic1,3)))
    Is = double(I(:,:,3,i));
    If1 = imfilter(Is,k1);
    If2 = imfilter(Is,k2);
    If3 = If1-If2*.5;
    Ic3(:,:,i) = If3;

    Is = double(I(:,:,2,i));
    If1 = imfilter(Is,k1);
    If2 = imfilter(Is,k2);
    If3 = If1-If2*.5;
    Ic2(:,:,i) = If3;


end

Ic2max = max(Ic2,[],3);
Ic3max = max(Ic3,[],3);

Ic2(Ic2<0) = 0;
Ic3(Ic3<0) = 0;

Irat = (Ic2-Ic3)./(Ic2+Ic3);
Irat(isnan(Irat)) = 0;

Iv = Irat.*Ic2;


%Irat = Irat*200+50;
IratMax = max(Irat,[],3);

for i = 1:size(Irat,3)
    image(Iv(:,:,i))
    drawnow
end


IvMax = max(Iv,[],3);
image(IvMax)


%% Threshold
thresh = 40;

It = Iv>=thresh;
image(sum(It,3)*8)

for i = 1:size(It,3)
    image(It(:,:,i)*1000)
    drawnow
end

Ilab = bwlabeln(It,6);
props = regionprops(Ilab,'BoundingBox','Area');
A = cat(1,props.Area);

big  = find(A>10000);

Icol1 = uint8(Ilab*0);
Icol2 = uint8(Ilab*0);
Icol3 = uint8(Ilab*0);

cmap = hsv(100);
cmap = cat(1,[ 0 0 0],cmap);
colIdx = ones(size(A)+1);
randCols = floor(rand(length(big),1)*99)+1;
colIdx(big+1) = randCols+1;

IcolIdx = colIdx(Ilab+1);
Icol1(:) = cmap(IcolIdx,1)*256;
Icol2(:) = cmap(IcolIdx,2)*256;
Icol3(:) = cmap(IcolIdx,3)*256;


IcSum1 = max(Icol1,[],3);
IcSum2 = max(Icol2,[],3);
IcSum3 = max(Icol3,[],3);

IcSum = cat(3,IcSum1,IcSum2,IcSum3);
image(uint8(IcSum))

sp1 = subplot(1,2,1);
image(sp1,IvMax)
sp2 = subplot(1,2,2);
image(sp2,IcSum);
colormap(sp1,gray(256))




%%




if 0
    %% Gamma
    Ig = I;
    for i = 1:size(I,4)
        for c = 3
            Is = Ig(:,:,c,i);
            Is2 = 256 .* (Is/256).^.6;
            image(Is2)

            Ig(:,:,c,i) = Is2;
        end
    end


    clf, colormap gray(256)
    for i = 1:size(I,4)

        Is = Ig(:,:,3,i);
        Is = (Is-20)* 1.4;
        Ig(:,:,3,i) = Is;

        image(Is)
        drawnow


    end
end

if 0

    clf, colormap gray(256)
    for i = 1:size(I,4)

        Is = Ig(:,:,:,i);
        fileName = sprintf('%s%s',TPN,nams{i});
        disp(fileName)
        imwrite(uint8(Is),fileName);


    end
end



