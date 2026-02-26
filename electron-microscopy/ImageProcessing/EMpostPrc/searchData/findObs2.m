colormap gray(255)
Istak = Istak;
Istak = Istak * 255/max(Istak(:));
maxIS = max(Istak(:));
for i = 1:size(Istak,3)
    i
    image(Istak(:,:,i)* 250/maxIS),pause
end

%%  threshold
% 
% for s = 1: size(Istak,3)
%    subplot(2,1,1)
%    image(Istak(:,:,s));
%    subplot(2,1,2)
%     
%     image((Istak(:,:,s)>50)*1000),pause
% end
% 

%% Search
pos = [];
for s = 1: size(Istak,3)
    s
    I = Istak(20 : end - 20,20 : end - 20,s);
    %I(I<0) = 0 ;
    subplot(2,2,1)
    image(I)


    [wI numWats] = bwlabel(I>40,4);
   
    peaks = I * 0;
    for i = 1:numWats
        ids = find(wI == i);
        vals = I(ids);
        getit = vals>=(max(vals)*.5);
        peaks(ids(getit))=vals(getit);
    end
    %image(peaks),pause(.01)

    [bwI nLab] = bwlabel(peaks);
    subplot(2,2,2)
    image(bwI*1000)

    stats = regionprops(bwI,I,'WeightedCentroid',...
        'Solidity','PixelIdxList','MajorAxisLength',...
        'MinorAxisLength','MeanIntensity');

    MinorAxisLength = [stats(:).MinorAxisLength];
    MajorAxisLength = [stats(:).MajorAxisLength];
    Solidity = [stats(:).Solidity];
    MeanIntensity = [stats(:).MeanIntensity];
    WC = [stats(:).WeightedCentroid];
    WC1 = WC(mod(1:length(WC),2)>0)';
    WC2 = WC(mod([1:length(WC)]+1,2)>0)';
    wCent = [WC1 WC2 ones(length(stats),1) * s ];
   
    
    centA =1;
    centSD = .05;
    cent = 1;
    gSolid = centA * exp(-.5 * ((Solidity-cent)/centSD).^2);
    %scatter(Solidity,gSolid)
    
    
    centA =1;
    centSD = 1;
    cent = 4;
    gMinAx = centA * exp(-.5 * ((MinorAxisLength-cent)/centSD).^2);
    scatter(MinorAxisLength,gMinAx)
    
    centA =1;
    centSD = 25;
    cent = 100;
    gMeanInt = centA * exp(-.5 * ((MeanIntensity-cent)/centSD).^2);
    gMeanInt(MeanIntensity>cent) = 1;
    scatter(MeanIntensity,gMeanInt)
    
    
    useStat = sum([2 * gMeanInt' 1 * gSolid' 1 * gMinAx'],2)/4 ;
    %useStat = gMeanInt .* gSolid .* gMinAx;
     
    showVal = bwI * 0;
    for i = 1: nLab
        ids = find(bwI == i);
        showVal(bwI == i) = useStat(i);
    end
    valStak(:,:,s) = showVal * 255;
    subplot(2,2,3)
    image(showVal * 255/max(showVal(:)))

    %hist(useStat)
    useReg = find(useStat>.75);
    
    showReg = showVal * 0;
    for i = 1: length(useReg); 
        ids = find(bwI == useReg(i));
        showReg(bwI == useReg(i)) =  showReg(bwI == useReg(i))+150;
    end
    subplot(2,2,4)
    image(showReg * 255/max(showReg(:))),pause(.1)
    
    pos = cat(1,pos, wCent(useReg,:));
    regStak(:,:,s) = showReg;
    colStak(:,:,1,s) = showReg;
    colStak(:,:,2,s) = I;
    colStak(:,:,3,s) = showVal;
end

subplot(1,1,1)
scatter3(pos(:,1),pos(:,2),pos(:,3));
%TPN = GetMyDir;
% %imwriteNp(TPN,Istak,'Istak')
imwriteNp(TPN,colStak,'colStak2')

mkdir([TPN '/pics/color'])
for s = 1:size(Istak,3)
   colP = uint8(colStak(:,:,1:3,s));
   imwrite(colP,[TPN 'pics/color/colStak' num2str(s) '.tif'],'compression','none')
    
end

