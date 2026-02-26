%%Find correlations between ROIs
%%What are correlations between nearby ROIs from different image planes

if 0
  
    global tis

    SPN = [glob.datDir 'Analysis\Data\preproc\'];
    load([SPN 'ptDat.mat']);
    load([SPN 'ROI.mat']);
    %load([SPN 'ROI2.mat']);
    load([SPN 'SOI.mat']);
    load([SPN 'RawOI.mat']);
    load([SPN 'NOI.mat']);
    load([SPN 'MOI.mat']);
    load([SPN 'COI.mat']);

end


%%
%sampRawBounds = [100 1500];
sampRawBounds = [100 300];

usedRois = [];
allL = [];
allR = [];
allE = [];
allry = [];
allrx = [];

for ry = 1:size(RawOI.useRois,1);
    for rx = 1:size(RawOI.useRois,2);

        pixResp = RawOI.PixResp{ry,rx};
        useRois = RawOI.useRois{ry,rx};
        rawPix = RawOI.rawPixResp{ry,rx};
        roiCid = SOI.cid(useRois);
        uRoi = 1:length(useRois) < [1:length(useRois)]';
        [roiY roiX]  = find(uRoi);

        plot(rawPix')

        pix = rawPix(:,sampRawBounds(1):sampRawBounds(2))';
        meanPix = mean(pix,1);
        difPix = pix-repmat(meanPix,[size(pix,1) 1]);
        R = corrcoef(pix);

        r2rLin = SOI.r2rLin(useRois,useRois);
        r2rEuc = SOI.r2rEuc(useRois,useRois);

        L = length(allR);
        N = sum(uRoi(:));
        allR(L+1:L+N) = R(uRoi);
        allL(L+1:L+N) = r2rLin(uRoi);
        allE(L+1:L+N) = r2rEuc(uRoi);
        allry(L+1:L+N) = ry;
        allrx(L+1:L+N) = rx;
        allRoiY(L+1:L+N) = useRois(roiY);
        allRoiX(L+1:L+N) = useRois(roiX);

    end
end

%% Display

sameCell = find(allL<inf);
difCell = find(allL==inf);

binW = 3;
hRange = [0:.1:150];


eSame = allE(sameCell);
rSame = allR(sameCell);
lSame = allL(sameCell);
eDif = allE(difCell);
lDif = allL(difCell);
rDif = allR(difCell);

clear reSame reDif rlDif rlSame
for i = 1:length(hRange)
    l = hRange(i);
    reSame(i) = mean(rSame((eSame>=(l-binW/2)) & (eSame<(l+binW/2))));
    rlSame(i) = mean(rSame((lSame>=(l-binW/2)) & (lSame<(l+binW/2))));
    reDif(i) = mean(rDif((eDif>=(l-binW/2)) & (eDif<(l+binW/2))));
    %rlDif(i) = mean(rDif((lDif>=(l-binW/2)) & (lDif<(l+binW/2))));
end

%%Make 3D surface relSame
binW2 = 5;
hRange2 =[0:1:150];

clear relSame
for i = 1:length(hRange2)
    l1 = hRange2(i);
    for o = 1:length(hRange2)
        l2 = hRange2(o);
        isE = (eSame>=(l1-binW2/2)) & (eSame<(l1+binW2/2));
        isL = (lSame>=(l2-binW2/2)) & (lSame<(l2+binW2/2));
        relSame(i,o) = mean(rSame(isE & isL));
        relSum(i,o) = sum(isE & isL);
    end
end
relSame(relSum<10) = nan;




clf
subplot(1,2,1)
hold on
scatter(eDif,rDif,2,'r','markerfacecolor','flat','markeredgealpha',.5,'MarkerFaceAlpha',.2)
scatter(lSame,rSame,2,'g','markeredgecolor',[0 .5 0],'markerfacecolor','flat','markeredgealpha',.5,'MarkerFaceAlpha',.2)
scatter(eSame,rSame,2,'b','markerfacecolor','flat','markeredgealpha',.5,'MarkerFaceAlpha',.2)
plot(hRange,reSame,'color',[0 0 0.4],'linewidth',2);
plot(hRange,reDif,'color',[.4 0 0],'linewidth',2);
plot(hRange,rlSame,'color',[0 0.4 0],'linewidth',2);

legend({' ',' ',' ','Euclidian distance between neurites of different cells',...
    'Euclidian distance between neurites of the same cell',...
    'Linear distance between neurites of the same cell'})
xlim([0 80])



subplot(1,2,2)
surf(hRange2,hRange2,relSame)
%daspect([1 1 .01])
daspect([1 1 100000])

view([90 -90])


if 0
    %fDir = uigetdir;
    TPN = uigetdir
    filename = [TPN '\autoCorr2_50-350'];
    set(gcf,'renderer','Painters')
    print('-depsc','-tiff','-r300', '-painters',[filename,'.eps'])
    
end


%% Find interesting comparisons
getDifNum = 100;
%%Get most correlated ROIs on different cells
rDif = allR(difCell);
[sortDif idx] = sort(rDif,'descend')
topDifs = difCell(idx(1:100));

topPairs = [allRoiY(topDifs)' allRoiX(topDifs)'];


SPN = [glob.datDir 'Analysis\Data\preproc\'];
ROIMask = load([SPN 'maskDat.mat']);
ROIMask = ROIMask.maskDat;

%%What is the shared influence for these correlated ROIs?

topPtIlab = ptDat(topPairs,1);
for i = 1:size(topPairs,1)
    
    clf

    lab = ptDat(topPairs(i,1),1);
    i1 = floor(lab/1000);
    i2 = mod(lab,1000);
    meanI1 = RawOI.meanIs{i1,i2} * .2;
    mask1 = (ROIMask(:,:,topPairs(i,1))>0)*30;
    colI1 = cat(3,meanI1,meanI1,meanI1);
    colI1(:,:,1) = colI1(:,:,1) + mask1;
    colI1(:,:,2) = colI1(:,:,2) - mask1;
    colI1(:,:,3) = colI1(:,:,3) - mask1;
    
    subplot(2,1,1)
    image(colI1)
    colormap gray(255)

    lab = ptDat(topPairs(i,2),1);
    i1 = floor(lab/1000);
    i2 = mod(lab,1000);
    meanI2 = RawOI.meanIs{i1,i2} * .2;
    mask2 = (ROIMask(:,:,topPairs(i,2))>0)*30;
    colI2 = cat(3,meanI2,meanI2,meanI2);
    colI2(:,:,1) = colI2(:,:,1) + mask2;
    colI2(:,:,2) = colI2(:,:,2) - mask2;
    colI2(:,:,3) = colI2(:,:,3) - mask2;

    subplot(2,1,2)
    image(colI2)
    colormap gray(255)
    pause

end









