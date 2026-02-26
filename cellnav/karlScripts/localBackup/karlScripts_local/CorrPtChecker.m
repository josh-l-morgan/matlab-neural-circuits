%% setup
analDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\';
datDir='Y:\karlsRetina\CellNavLibrary_IxQ\Analysis\Data\preproc\';
physDir='W:\CorrectRegistratedFiles\CorrectRegistratedFiles\';
ptsFileName='ptDat.mat';

initLoadBool=1;
if initLoadBool
    global regionNames
    global pts
    global dsObj
    global tis
    global corrPts
    global physImgs
    global vgcMat
    load([datDir ptsFileName]);
    load([analDir 'Merge\dsObj.mat']);
    load([analDir 'Analysis\tis.mat']);
    regionNames=[1005:1010 2001:2006];
    %% load
    [corrPts,physImgs,vgcMat]=initLoad(2,physDir,regionNames,pts,dsObj,tis);
end
ptDat=sortrows(ptDat,[1,4,5]);
pts.ptDat=ptDat;
pts.ptDat=sortrows(pts.ptDat,[1,4,5]);
corrPts=refreshPts(pts,regionNames);
allTFs=refreshTfs(corrPts);

%% make a pretty figure
i=9;
%figA=plotPhys(corrPts,regionNames,i,physImgs);
%figB=plotPlane(corrPts,i);
%figC=projectEMontoPhys(corrPts,regionNames,i,physImgs,vgcMat);
newPts=addPhysPts(pts,i,regionNames,physImgs,vgcMat,figA,figB,figC);

%% We've got the func
function newPts = addPhysPts(pts,i,regionNames,physImgs,vgcMat,figA,figB,figC)
corrPts=refreshPts(pts,regionNames);
allTFs=refreshTfs(corrPts);
curPts=corrPts{i};
newPts=pts;
figA=plotPhys(corrPts,regionNames,i,physImgs);
%figure(figName);
hold on
while true
    pts=newPts;
    figure(figA);
    hold on
    [x,y,butt]=ginput();
    if isempty(x)
        break
    elseif butt==3
        'removing previous pt'
        %remove pt
    else
        scatter(x,y,50,'r+');
        [vastXds,vastYds]=allTFs{i}.transformPointsForward(x/10,y/10);
        vastimate=[vastXds,vastYds]*25;
        curPlanePts=curPts(:,6:8);
        Bp=[curPlanePts(:,1), curPlanePts(:,2), ones(size(curPlanePts,1),1)] \ curPlanePts(:,3);
        vastZ=Bp(1)*vastimate(1) + Bp(2)*vastimate(2) + Bp(3)*ones(size(vastimate(1)));
        vastimate = [vastimate(1),vastimate(2),vastZ];
        vastOut=round(vastimate);
        clipboard('copy',vastOut);
        newVast=input('newVast?');
        targetCid=input('targetVGC?');        
        curNewCorr(1)=regionNames(i);
        curNewCorr(3)=targetCid;
        curNewCorr(4:5)=[round(x/10)-10 round(y/10)-10];
        curNewCorr(6:8)=newVast;
        %curNewCorr(q,8)=round(sum(emImageOffsets{i}([3 6]))/2*2.5);
        pts.ptDat=vertcat(pts.ptDat,curNewCorr);
        %pts.ptDat(
        corrPts=refreshPts(pts,regionNames);
        allTFs=refreshTfs(corrPts);
        %update planeFig
        %update physFig
        %figure(figName);
        newPts=pts;
        %figName=plotPhys(corrPts,regionNames,i,physImgs);
        figB=plotPlane(corrPts,i);
        figC=projectEMontoPhys(corrPts,regionNames,i,physImgs,vgcMat);
    end
end
end

function outFig = projectEMontoPhys(corrPts,regionNames,i,physImgs,vgcMat)
outFig=plotPhys(corrPts,regionNames,i,physImgs);
%vgcTest=sum(vgcMat,3);
allTFs=refreshTfs(corrPts);
curPts=corrPts{i};
zrange=[min(curPts(:,8)) max(curPts(:,8))]/2.5;
[osX,osY]=allTFs{i}.transformPointsForward([0,256,256,0],[0,0,32,32]);
debug=0;
if debug
    figure(); imshow(vgcTest);
    hold on
    title('Project PhysRegion onto native EM');
    plot(osX([1 2 3 4 1]),osY([1 2 3 4 1]),'r-');
end
[iosX,iosY]=allTFs{i}.transformPointsForward([-10,266],[-10,42]);
%cornY=min(osY);
%cornX=min(osX);
%offSet=[osX,osY]*2.5;
invertTF=invert(allTFs{i});
scaleTF=affine2d;
scaleTF.T(1,1)=10;scaleTF.T(2,2)=10;
invertUPTF=invertTF;
invertUPTF.T=invertTF.T*scaleTF.T;
%vgcSubMat=sum(vgcMat(iosY(1):iosY(2),iosX(1):iosX(2),zrange(1)-25:zrange(2)+25),3);
vgcSubMat=sum(vgcMat(:,:,zrange(1)-25:zrange(2)+25),3);
%figure(); imshow(vgcSubMat);
[testImg,testRef]=imwarp(vgcSubMat,invertUPTF);
%[emOrigX,emOrigY]=invertUPTF.transformPointsForward(osX(1),osY(1));
[emOrigX,emOrigY]=testRef.worldToIntrinsic(-100,-100);
%figure(); imshow(testImg);
crpImg=testImg(floor(emOrigY(1)):floor(emOrigY(1))+519,floor(emOrigX(1)):floor(emOrigX(1))+2759);
%figure(); imshow(crpImg);
%emImgRef=imref2d([2560,320],[100,2660],[100,420]);
physImgUp=imresize(physImgs{i},10)/255;
%figure(); imshow(physImgUp);
regdImg=zeros(520,2760,3);
regdImg(:,:,1)=crpImg;
regdImg(100:419,100:2659,2)=physImgUp;
hold on
imshow(regdImg);
end

function tfs = refreshTfs(corrPts)
tfs={};
for i=1:length(corrPts)
    curPts=corrPts{i};
    %fit the plane
    vastCoords=curPts(:,6:8);
    vastCoordsDS=vastCoords./[25 25 2.5];
    physPts=horzcat(curPts(:,4:5),repmat(1,[size(curPts,1),1]));
    tf3d=fitgeotrans(physPts(:,[1 2]),vastCoordsDS(:,[1 2]),'affine');
    tfs{i}=tf3d;
end
end

function corrPts = refreshPts(pts,regionNames)
corrPts=cell(12,1);
for w=1:length(regionNames)
    curPts=pts.ptDat(pts.ptDat(:,1)==regionNames(w),:);
    corrPts{w}=curPts;
end
end

function outFig=plotPlane(corrPts,i)
curPlanePts=corrPts{i}(:,6:8);
Bp=[curPlanePts(:,1), curPlanePts(:,2), ones(size(curPlanePts,1),1)] \ curPlanePts(:,3);
[X,Y]=meshgrid(linspace(25000,45000,10),linspace(30000,35000,16));
%zest=(-double(x)*GCLplane.Parameters(2)-double(y)*GCLplane.Parameters(3)-GCLplane.Parameters(4))/GCLplane.Parameters(1);
Z = Bp(1)*X + Bp(2)*Y + Bp(3)*ones(size(X));
estZ=Bp(1)*curPlanePts(:,1) + Bp(2)*curPlanePts(:,2) + Bp(3)*ones(size(curPlanePts(:,1)));
linErr=abs(estZ-curPlanePts(:,3));
outFig=figure();
scatter3(curPlanePts(:,1),curPlanePts(:,2),curPlanePts(:,3));
hold on
for w=1:length(curPlanePts(:,1))
    plot3([curPlanePts(w,1) curPlanePts(w,1)], ...
        [curPlanePts(w,2) curPlanePts(w,2)], ...
        [curPlanePts(w,3) estZ(w)],'r-');
end
end

function outFig=plotPhys(corrPts,regionNames,i,physImgs)
curPts=corrPts{i};
imgBord=100;
outFig=figure();

physImgUp=imresize(physImgs{i},1);
physImgRef=imref2d(size(physImgUp),[100,2660],[100,420]);
emptyBG=zeros(520,2760);
imshow(emptyBG);
hold on
title(string(regionNames(i)));
imshow(physImgUp,physImgRef,[0 255]);
scatter(curPts(:,4)*10+imgBord,curPts(:,5)*10+imgBord,25,'ro','filled');
end

function [corrPts,physImgs,vgcMat]=initLoad(level,physDir,regionNames,pts,dsObj,tis)
%level 0 = no load
%level 1 = load vgcVox
%level 2 = load vgcVox + Phys
if level>0
corrPts=cell(12,1);
if level>1
physImgs=cell(12,1);
end
for i=1:length(regionNames)
    if level>1
    curPhysImgFilename=strcat(['Ai148_129SVG3_Translation_122618_' num2str(regionNames(i)) '.mat']);
    curPhysImgDat=load([physDir curPhysImgFilename]);
    curPhysMean=mean(curPhysImgDat.I,3);
    curPhysMeanImg=curPhysMean*(256/prctile(curPhysMean(:),95));
    physImgs{i}=curPhysMeanImg;
    curPhysImFilename=['Physimage_' num2str(regionNames(i)) '.png'];
    %imwrite(curPhysMeanImg/255,['Y:\MATLAB\cellNav\cellNav_0.70\karl\images\emPhys\' curPhysImFilename]);
    end
    curPts=pts.ptDat(pts.ptDat(:,1)==regionNames(i),:);
    corrPts{i}=curPts;
end
vgcCidList=[2 3 4 5 10 11 13 14 20 ];
cidVoxs=getCidVox(vgcCidList,1,dsObj,tis);
emptyMat=zeros(3200,3200,500,'uint8');
vgcMat=emptyMat;
for i=1:length(cidVoxs)
    curVoxList=cidVoxs{i};
    curVoxList=curVoxList(sum(curVoxList,2)>2,:);
    vgcMat(sub2ind(size(vgcMat),curVoxList(:,1),curVoxList(:,2),curVoxList(:,3)))=i;
end
end
end