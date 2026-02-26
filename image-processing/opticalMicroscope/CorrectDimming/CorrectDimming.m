

SPN = 'H:\LxA\mx\TifStack\'
SPN = 'H:\LxA\mx\redGreenGray\'
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

if 0

    clf, colormap gray(256)
    for i = 1:size(I,4)

        Is = Ig(:,:,:,i);
        fileName = sprintf('%s%s',TPN,nams{i});
        disp(fileName)
        imwrite(uint8(Is),fileName);


    end
end



