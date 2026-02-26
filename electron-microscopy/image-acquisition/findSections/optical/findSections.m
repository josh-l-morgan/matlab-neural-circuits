clear all
tic
colormap gray(255)
[TFN TPN] = GetMyFile;


%%
I = imread([TPN TFN]);
image(I)
%I = I(100:600,300:400,:);
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
canI = edge(I,'Canny',[.01 .011],1);
[canI cNum] = bwlabel(canI,8);
image(fitH(canI)),colormap colorcube(255)
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
colI(:,:,3) = canI*100;
image(colI),pause(1)

[lI oNum] = bwlabel(oI>0,8); %%label objects
%%

eoI = I*0;
image(fitH(eoI>1))
for o = 1: oNum
    [oy ox] = find(lI == o);
    sI = I(min(oy):max(oy),min(ox):max(ox));
    idO = sub2ind(size(sI),oy-min(oy)+1,ox-min(ox)+1);
    esI = edge(sI,'Canny',[.01 .5],1);
    lE = bwlabel(esI,8);
    goodEdges = lE(idO);
    goodEdges = unique(goodEdges(goodEdges>0));
    badEdges = setdiff(1:max(lE(:)),goodEdges);
    for e = 1:length(badEdges)
        lE(lE==badEdges(e))=0;
    end
    eoI(min(oy):max(oy),min(ox):max(ox))= eoI(min(oy):max(oy),min(ox):max(ox))+ (lE>0)*o;
    
end
image(eoI * 100),pause(.1)
%% combine objects with edges

eI = canI * 0;
for o = 1:oNum
    edges = canI(lI==o);
    edges = unique(edges(edges>0))
    for e = 1:length(edges)
        eI(canI == edges(e)) = o;
    end
    
end
image(eI),colormap colorcube(255)

eIa = eI * 0;
eProps = regionprops(eI,'Orientation','Centroid','ConvexHull');
for i = 1:length(eProps)
    cPoly = eProps(i).ConvexHull;
   eIa(poly2mask(cPoly(:,1),cPoly(:,2),ys,xs)) = 1;
end
image(eIa*1000),pause

%%

rProps = regionprops(oI,cI2,'Area','EulerNumber','Orientation',...
    'WeightedCentroid','meanIntensity','majorAxisLength',...
    'MinorAxisLength','Perimeter','Extent','Extrema','ConvexHull',...
    'PixelIdxList','BoundingBox');
% toc
% 
% showI = I * 0;
% for i = 1:length(rProps)
%     showI(sub2ind(size(I),rProps(i).Extrema(:,1),rProps(i).Extrema(:,2)))=100;
% end
% colI(:,:,1) = showI;
% image(colI)
% 
% rProps(3).Extrema

%% Look for verts
gKern = zeros(100:100);
gKern(:,50) = 100;
vI = fastCon(I,gKern);
image(fitH(vI))

%%  Shrink to edges
for i = 1:max(oI(:))
    getI = canI*100+(oI==i)*50;
    
   bBox = rProps(i).BoundingBox
   bBox(3) = fix(bBox(1)+ bBox(3));
   bBox(4) = fix(bBox(2) + bBox(4));
   bBox = fix(trim(bBox,size(I')));
%    scanI = canI(bBox(2):bBox(4),bBox(1):bBox(3));
%    grabO = oI(bBox(2):bBox(4),bBox(1):bBox(3));
%    getI(bBox(2):bBox(4),bBox(1):bBox(3)) = getI(bBox(2):bBox(4),bBox(1):bBox(3)) + scanI*1000;
%   
   
   %getI(oI==i) = I(oI==i);
  % scanI = edge(getI,'Canny',[],1);
   %image(getI + scanI*1000);
   
   image(fitH(getI))
 pause

end

%% Class

scatter(minAx,majAx,'.')
