clear all
tic
colormap gray(255)
[TFN TPN] = GetMyFile;


%%
I = imread([TPN TFN]);
image(I)
%I = I(1:400,1:400,:);
image(I(:,:,1)./I(:,:,2))
[ys xs] = size(I);

I = mean(I,3);
I = max(I(:)) - I;
colormap gray(255)
image(I)

gKern = gaus3d([5 5 1],3,[1 1 1]);
cI = fastCon(I,gKern);
dI = (I - cI);
image(dI*50+100)
image(cI)

%%
%canI = edge(I,'Canny',[.08 .1],1);
vI = edge(I,'Canny',[.03 .1],2);
cvI = imclose(vI,strel('disk',10));
image(fitH(vI))
image(fitH(cvI))
[cvI cNum] = bwlabel(cvI==0,4);
vProps = regionprops(cvI,'Orientation','MajorAxisLength',...
    'MinorAxisLength','ConvexHull');
clear minAx maxAx ori
for i = 1:length(vProps)
    minAx(i) = vProps(i).MinorAxisLength
    maxAx(i) = vProps(i).MajorAxisLength
    ori(i) = vProps(i).Orientation;
end
scatter(minAx,maxAx,'.')

tape = find(~((minAx>80) & (minAx<120) & ((maxAx./minAx)>4)));
for i = 1:length(tape)
    cvI(cvI==tape(i)) =0;
end
image(fitH(cvI))

%%
gKern = gaus3d([10 10 1],5,[1 1 1]);
cI2 = fastCon(I,gKern);
image(cI2)

dI2 = max(cI2(:))-cI2;
%dI2 = dI2 * 100/median(dI2(dI2>0));
%dI2(~I) = -Inf;%min(dI2(:));
%gKern = gaus3d([15 15 5],2,[1 1 3]);
%cI = fastCon(dI2,gKern);
dI2 = imhmin(dI2,2); %
%dI2(dI2<-10) = -10;

wI = watershed(dI2,8);
image(wI)
image((wI==0)*1000 + dI2);

%%
oI = I * 0;
for i = 1:max(wI(:))
    tI =  cI2;
    tI(wI~= i) = 0;
    maxT = max(tI(:));
    minT = min(tI(tI>0));
    pList = find(bwperim(tI));
    seed = find(tI == maxT,1);
    %image(tI),pause(.01)
    %hist(dI2(wI==i),1:1:200)
    for t = maxT:-1:minT
        [lI nl] = bwlabel(tI>=t);
        if  ~isempty(find(lI(pList)==lI(seed)))
            break
        end
        oI(lI==lI(seed)) = oI(lI==lI(seed))+1;
    end
end
image((oI>0)*1000)
colI = uint8(cat(3,oI,oI,oI) * 50);
colI(:,:,1) = dI * 50 +100;
image(colI),pause(1)

oI(cvI == 0) = 0;
[lI oNum] = bwlabel(oI>0,8); %%label objects
%%

eI = I*0;
image(fitH(eI>1))
for o = 1: oNum
    [oy ox] = find(lI == o);
    sI = I(min(oy):max(oy),min(ox):max(ox));
    idO = sub2ind(size(sI),oy-min(oy)+1,ox-min(ox)+1);
    esI = edge(sI,'Canny',[.001 .2],1.5);
    lE = bwlabel(esI,8);
    goodEdges = lE(idO);
    goodEdges = unique(goodEdges(goodEdges>0));
    badEdges = setdiff(1:max(lE(:)),goodEdges);
    for e = 1:length(badEdges)
        lE(lE==badEdges(e))=0;
    end
    eI(min(oy):max(oy),min(ox):max(ox))= eI(min(oy):max(oy),min(ox):max(ox))+ (lE>0)*o;
    
end
image(eI),colormap colorcube(255),pause(.1)
%% combine objects with edges

eIa = eI * 0;
eProps = regionprops(eI,'Orientation','Centroid','ConvexHull');
for i = 1:length(eProps)
    cPoly = eProps(i).ConvexHull;
   eIa(poly2mask(cPoly(:,1),cPoly(:,2),ys,xs)) = 1;
end
image(eIa*1000),pause(.1)
colI(:,:,3) = fitH(eIa);
colI(:,:,2) = fitH(I);
colI(:,:,1) = fitH(cvI>0);
image(colI),pause(.1)

[lI oNum] = bwlabel(eIa>0,4); %%label objects
%% select sections

rProps = regionprops(lI,cI2,'Area','EulerNumber','Orientation',...
    'WeightedCentroid','MeanIntensity','MajorAxisLength',...
    'MinorAxisLength','Perimeter','Extent','Extrema','ConvexHull',...
    'PixelIdxList','BoundingBox');

minAx = [rProps.MinorAxisLength];
majAx = [rProps.MajorAxisLength];
scatter(minAx,majAx,'.')

T = cluster([minAx;majAx],'maxclust',2)
T = kmeans([minAx;majAx],2);
[IDX,C,sumd,D] = kmeans([minAx]',2)

scatter(minAx(IDX == 1),majAx(IDX == 1),'.','r');
hold on
scatter(minAx(IDX == 2), majAx(IDX == 2),'.','b');
hold off

useK = find(IDX == find(C == min(C)));
for i = 1: length(useK)
   lI(lI == useK(i)) = 0; 
end

colI(:,:,3) = fitH(lI>0);
image(colI),pause(.1)

[lI oNum] = bwlabel(lI>0,4);
rProps = regionprops(lI,'Centroid','ConvexHull');



%%