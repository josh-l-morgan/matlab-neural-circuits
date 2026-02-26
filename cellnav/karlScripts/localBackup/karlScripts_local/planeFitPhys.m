%% load data
dataDir='Y:\MATLAB\cellNav\cellNav_0.67\karl\';
ptDat=load([dataDir 'corrExport.csv']);
roiList=unique(ptDat(:,1));
imMult=3;
markSize=10;
linWidth=2;
colmap=[1 0 0 ; 0 0 1 ; 0 1 0 ; 1 1 0 ; 1 0 1 ; 0 1 1];
corrPtError={};
deplanes=struct();

%% establish that affine is good for these things.

for i=length(roiList):-1:1
    curPts=ptDat(ptDat(:,1)==roiList(i),:);
    curPhysPts=curPts(:,[4 5]);
    curEMPts=curPts(:,[6 7]);
    tform=fitgeotrans(curEMPts, curPhysPts, 'affine');
    deplanes(i).aff2dTform=tform;
    EMtfd=transformPointsForward(tform,curEMPts);
    dists=sqrt(sum(((EMtfd-curPhysPts).^2),2));
    deplanes(i).aff2dError=dists;
end
%tform = fitgeotrans(moving, fixed, type);

%% fit de planes

figure();
hold on
B=zeros(3,1);
for i=1:length(roiList)
    P=ptDat(ptDat(:,1)==roiList(i),[8 6 7]);
    B(:) = [P(:,3), P(:,2), ones(size(P,1),1)] \ P(:,1);
    %curPlane=struct();
    %curPlane.Parameters=[-1 B(2,1) B(1,1) B(3,1)];
    deplanes(i).avgZ=mean(P(:,1));
    deplanes(i).points=P;
    deplanes(i).B=B;
    scatter3(P(:,2),P(:,3),P(:,1), 25, colmap(mod(i,6)+1,:));
    if i<7
        y = [27000:1000:30000];
    else
        y = [32000:1000:35000];
    end
    x = [25000:1000:45000];
    [X,Y]=meshgrid(x,y);
    z = (B(2)*X + B(1)*Y + B(3));
    surf(x, y, z);
    %roiList(i)
    planeZ=(B(2)*P(:,2) + B(1)*P(:,3) + B(3));
    ZError=P(:,1)-planeZ;
    absZError=abs(ZError);
    deplanes(i).Zerror=ZError;
    %i=i+1;
end
title('Correspondence Pts and Fitted Planes');

%% get the average plane and find the z averages

goodPlanes=[1:6];
totB=zeros(3,1);
denom=0;
for i=1:length(goodPlanes)
    curPlane=goodPlanes(i);
    totB=totB+length(deplanes(curPlane).Zerror)*deplanes(curPlane).B;
    denom=denom+length(deplanes(curPlane).Zerror);
end
avgB=totB/length(goodPlanes)/denom;
avgB(3)=600;
x = [25000:1000:45000];
y = [25000:1000:29000];
[X,Y]=meshgrid(x,y);
z = (avgB(2)*X + avgB(1)*Y + avgB(3));
surf(x, y, z);
avgB

goodPlanes=[7:12];
totB=zeros(3,1);
denom=0;
for i=1:length(goodPlanes)
    curPlane=goodPlanes(i);
    totB=totB+length(deplanes(curPlane).Zerror)*deplanes(curPlane).B;
    denom=denom+length(deplanes(curPlane).Zerror);
end
avgB=totB/length(goodPlanes)/denom;
avgB(3)=600;
x = [25000:1000:45000];
y = [31000:1000:35000];
[X,Y]=meshgrid(x,y);
z = (avgB(2)*X + avgB(1)*Y + avgB(3));
surf(x, y, z);
avgB


%% asdf
figure();
hold on
title('Correspondence Pts on Summed Functional Images');
plotLoc=[1 3 5 7 9 11 2 4 6 8 10 12];
for i=1:length(roiList)
    curROIstring=num2str(roiList(i));
    %curTiff=Tiff([dataDir 'images\AVG_Ai148_129SVG3_122618_' curROIstring '-1.tif'], 'r');
    imDat=imread([dataDir 'images\SUM_Ai148_129SVG3_122618_' curROIstring '-1.tif']);
    deplanes(i).imDat=imDat;
    %imDat=read(curTiff);
    imDatX=imDat*imMult;
    subplot(6,2,plotLoc(i));
    imshow(imDatX);
    hold on;
    title(curROIstring);
    curPts=ptDat(ptDat(:,1)==roiList(i),:);
    for j=1:length(curPts(:,1))
        plot(curPts(j,4),curPts(j,5),'w+','MarkerSize', markSize, 'LineWidth', linWidth);
        plot(curPts(j,4),curPts(j,5),'mo','MarkerSize', deplanes(i).aff2dError(j)*2, 'LineWidth', linWidth);
        plot(curPts(j,4),curPts(j,5),'co','MarkerSize', abs(deplanes(i).Zerror(j)+1/50), 'LineWidth', linWidth);
    end
end
hold off

ha=get(gcf,'children');
for k=1:length(ha)
    set(ha(k),'position',[floor((k-1)/6)/2 mod((k-1),6)/6+.01 .49 .1])
end
%set(ha(7),'position',[.5 .17 .49 .1])
%set(ha(3),'position',[.5 .5 .4 .4])
%set(ha(4),'position',[.1 .5 .4 .4])


%% try the first figure again with average plane slopes
% at average Z displacement

goodPlanes=[1:12];
totB=zeros(3,1);
denom=0;
for i=1:length(goodPlanes)
    curPlane=goodPlanes(i);
    totB=totB+length(deplanes(curPlane).Zerror)*deplanes(curPlane).B;
    denom=denom+length(deplanes(curPlane).Zerror);
end
tavgB=totB/length(goodPlanes)/denom;

%% avg planes plot

figure();
hold on
B=zeros(3,1);
for i=1:length(roiList)
    P=ptDat(ptDat(:,1)==roiList(i),[8 6 7]);
    scatter3(P(:,2),P(:,3),P(:,1), 25, colmap(mod(i,6)+1,:));
    if i<7
        y = [27000:1000:30000];
    else
        y = [32000:1000:35000];
    end
    x = [25000:1000:45000];
    [X,Y]=meshgrid(x,y);
    B=tavgB;
    B(3)=deplanes(i).avgZ;
    z = (B(2)*X + B(1)*Y + B(3));
    surf(x, y, z);
    %roiList(i)
    planeZ=(B(2)*P(:,2) + B(1)*P(:,3) + B(3));
    ZError=P(:,1)-planeZ;
    absZError=abs(ZError);
    deplanes(i).avgZerror=ZError;
    %i=i+1;
end
title('Correspondence Pts and avg planes');


%% asdf
figure();
hold on
for i=1:length(roiList)
    P=ptDat(ptDat(:,1)==roiList(i),[8 6 7]);
    scatter(repmat(i,length(P(:,1))),P(:,1), 25, colmap(mod(i,6)+1,:));
end
xticks([1:12]);
xticklabels(cellstr(num2str(roiList)));
title('Z distribution of Ephys Correspondence Points by Imaging Field');

figure();
hold on
for i=1:length(roiList)
    P=ptDat(ptDat(:,1)==roiList(i),[8 6 7]);
    scatter(repmat(i,length(P(:,1))),deplanes(i).Zerror, 25, colmap(mod(i,6)+1,:));
end
xticks([1:12]);
xticklabels(cellstr(num2str(roiList)));
title('Z error of Correspondence Points');

figure();
hold on
for i=1:length(roiList)
    P=ptDat(ptDat(:,1)==roiList(i),[8 6 7]);
    scatter(repmat(i,length(P(:,1))),deplanes(i).avgZerror, 25, colmap(mod(i,6)+1,:));
end
xticks([1:12]);
xticklabels(cellstr(num2str(roiList)));
title('Dist to Avg Planes of Correspondence Points');

%% Get bounding boxes for each ROI
boundBox=[1 1 ; 256 1 ; 1 32 ; 256 32];
for i=1:length(roiList)
    B=tavgB;
    B(3)=deplanes(i).avgZ;
    BBtfd=transformPointsInverse(deplanes(i).aff2dTform,boundBox);
    z = (B(2)*boundBox(:,1) + B(1)*boundBox(:,2) + B(3));
    deplanes(i).bbox=[BBtfd z];
    
end


%% test
if 0
    P=Locs;
    B(:,1) = [P(:,3), P(:,2), ones(size(P,1),1)] \ P(:,1);
    curPlane=struct();
    curPlane.Parameters=[-1 B(2,1) B(1,1) B(3,1)];
    testFig=1;
    if testFig==1
        figure();
        hold on
        scatter3(Locs(:,2),Locs(:,3),Locs(:,1),'bo');
        %scatter3(INLbord(1:10:end,2),INLbord(1:10:end,3),INLbord(1:10:end,1),'ro');
        for i=1
            x = [20000:1000:40000];
            y = [20000:1000:40000];
            [X,Y]=meshgrid(x,y);
            z = (B(2,i)*X + B(1,i)*Y + B(3,i));
            surf(x, y, z);
        end
        %scatter3(verts(:,2),verts(:,3),verts(:,1));
        
        testResults=zeros(10,4);
        if 0
            for i=1:10
                curPlane=pcfitplane(GCLptCloud,1,'MaxNumTrials',50000);
                testResults(i,:)=curPlane.Parameters(:);
                x = [50:10:200];
                y = [50:10:200];
                [X,Y]=meshgrid(x,y);
                curParams=curPlane.Parameters(:);
                if curParams(4)<0
                    curParams=curParams*-1;
                end
                z = (curParams(2)*X + curParams(3)*Y + curParams(4));
                surf(x, y, z);
            end
        end
        
        if 0
            for i=1:10
                INLplane=pcfitplane(INLptCloud,1,'MaxNumTrials',50000);
                testResults(i,:)=INLplane.Parameters(:);
                x = [50:10:200];
                y = [50:10:200];
                [X,Y]=meshgrid(x,y);
                curParams=INLplane.Parameters(:);
                if curParams(4)<0
                    curParams=curParams*-1;
                end
                z = (curParams(2)*X + curParams(3)*Y + curParams(4));
                surf(x, y, z);
            end
        end
        
    end
end