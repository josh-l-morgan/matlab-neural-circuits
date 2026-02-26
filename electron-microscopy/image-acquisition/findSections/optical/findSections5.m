clear all

colormap gray(255)
[TFN TPN] = GetMyFile;
tic
estimateWidth = 25;
eW = estimateWidth;


%%  Get Image
I = imread([TPN TFN]);
I = I(200:600,200:500,:);
[ys xs] = size(I);

I = mean(I,3);
I = max(I(:)) - I;

subplot(1,2,1)
image(fitH(I)),pause(.01)
subplot(1,2,2)

gKern = gaus3d([10 10 1],3,[1 1 1]);
smI = fastCon(I,gKern);
gKern = gaus3d([64 64 1],40,[1 1 1]);
bgI = fastCon(I,gKern);
difI = (smI - bgI);
image(fitH(difI))

%% Find Tape
vI = edge(difI,'Canny',[.001 .08],2);
image(fitH(vI))
cvI = imclose(vI,strel('disk',10));
image(fitH(cvI))
[cvI cNum] = bwlabel(cvI==0,4);
vProps = regionprops(cvI,'Orientation','MajorAxisLength',...
    'MinorAxisLength','ConvexHull');
clear minAx maxAx ori

minAx = [vProps.MinorAxisLength];
maxAx = [vProps.MajorAxisLength];
ori = [vProps.Orientation];
scatter(minAx,maxAx,'.')

tape = find(((minAx>50) & (minAx<110) & ((maxAx./minAx)>4)));
tpI = I * 0;
for i = 1:length(tape)
      cPoly = vProps(tape(i)).ConvexHull;
    tpI(poly2mask(cPoly(:,1),cPoly(:,2),ys,xs)) = 1;
end
image(fitH(tpI)),pause(.1)

%% Find sections

xI = edge(difI,'Canny',[.0001 .05],2);
image(fitH(xI))
%%

gKern = gaus3d([10 10 1],5,[1 1 1]);
cI2 = fastCon(I,gKern);

dI2 = cI2;
dI2 = max(dI2(:))-dI2;
dI2 = imhmin(dI2,3);

wI = watershed(dI2,8);
image((wI==0)*1000 + dI2);

%Arbitrary trimming
image((wI==0) * 100 + (xI)*100)
xI(wI==0) = 0;



%% filter edges
xI(~tpI) = 0;
[lI oNum] = bwlabel(xI,8);
image(lI * 1000)

xProps = regionprops(lI,'MinorAxisLength','MajorAxisLength',...
    'PixelIdxList','ConvexHull');
minAx = [xProps.MinorAxisLength];
majAx = [xProps.MajorAxisLength];
bad = find((minAx < 5) | (majAx > 100));

for b = 1:length(bad)
   lI(lI==bad(b))=0; 
end
%remove small breaks
cvI = imclose(lI>0,strel('disk',4));
image(fitH(cvI))

lI = bwlabel(lI,8);


xProps = regionprops(lI,'MinorAxisLength','MajorAxisLength',...
    'PixelIdxList','ConvexHull');
valI = xI *0;
good = 1:length(xProps);
for x = 1:length(good)
    cPoly = xProps(good(x)).ConvexHull;
    idPoly = find(poly2mask(cPoly(:,1),cPoly(:,2),ys,xs));
    valI(idPoly) = mean(difI(idPoly));
end
image(fitH(valI))



lI = bwlabel(lI,8);
xProps = regionprops(lI,'MinorAxisLength','MajorAxisLength',...
    'PixelIdxList')

%% 
% xI = edge(I,'Canny',[.00001 .02],5);
% image(fitH(xI))
% cxI = imclose(xI,strel('disk',5));
% image(fitH(cxI))
% [cxI cNum] = bwlabel(cxI==0,4);
% xProps = regionprops(cxI,'Orientation','MajorAxisLength',...
%     'MinorAxisLength','ConvexHull');
% clear minAx maxAx ori
% 
% minAx = [vProps.MinorAxisLength];
% maxAx = [vProps.MajorAxisLength];
% ori = [vProps.Orientation];
% scatter(minAx,maxAx,'.')
% 
% tape = find(((minAx>80) & (minAx<110) & ((maxAx./minAx)>4)));
% tpI = I * 0;
% for i = 1:length(tape)
%       cPoly = vProps(tape(i)).ConvexHull;
%     tpI(poly2mask(cPoly(:,1),cPoly(:,2),ys,xs)) = 1;
% end
% image(fitH(tpI)),pause(.1)
% 
% 
%  % Watershed
% gKern = gaus3d([10 10 1],5,[1 1 1]);
% cI2 = fastCon(I,gKern);
% 
% dI2 = max(cI2(:))-cI2;
% dI2 = imhmin(dI2,1);
% 
% wI = watershed(dI2,4);
% image((wI==0)*1000 + dI2);
% 
% % %% find wells
% profile on
% oI = I * 0;
% outside= find(bwperim(I+1));
% for i = 1:max(wI(:))
%     if ~mod(i,100)
%         sprintf('%d of %d',i,max(wI(:)))
%     end
%     tI =  cI2;
%     tI(wI~= i) = 0;
%     maxT = max(tI(:));
%     minT = min(tI(tI>0));
%     pList = find(bwperim(tI));
%     %pList = setdiff(pList,outside);
%     %tI(pList) = cI(pList);
%     seed = find(tI == maxT,1);
% %     for t = maxT:-1:minT
% %         [lI nl] = bwlabel(tI>=t);
% %          oI(lI==lI(seed)) = oI(lI==lI(seed))+1;
% %         if  ~isempty(find(lI(pList)==lI(seed)))
% %             break
% %         end
% %        
% %     end
% %     
%     minT = median(cI(pList));
%     thresh = maxT - (maxT-minT)/2;
%     [lT tNum] = bwlabel(tI>thresh);
%     %image(fitH(lT)),pause
%     oI(lT==lT(seed)) = tI(lT==lT(seed))-thresh;
%   
% end
% oI(oI<0)=0;
% % 
% %   image(fitH(oI)),pause
% % %%Trim objects
% % [lI oNum] = bwlabel(oI>0,4); %%label objects
% % for o = 1:oNum
% %    idx = find(lI==o); 
% %    idVal = oI(idx);
% %    idVal = idVal-max(idVal(:)) * .5;
% %    idVal(idVal<0) = 0;
% %    oI(idx) = idVal;
% % end
% colI = uint8(cat(3,oI,oI,oI) * 50);
% colI(:,:,1) = dI * 50 +100;
% image(colI),pause(1)
% 
% oI(tpI == 0) = 0;
% image(fitH(oI))
% profile off
% %% Describe objects.
% 
% 
% [lI oNum] = bwlabel(oI>0,8); %%label objects
% oProps = regionprops(lI,oI,'MinorAxisLength','MajorAxisLength');
% % [center,U,objFcn] = fcm(oDat,2);
% minAx = [oProps.MinorAxisLength]';
% tooSmall = find(minAx < eW/3);
% for i = 1:length(tooSmall)
%    oI(lI==tooSmall(i)) = 0; 
% end
% [lI oNum] = bwlabel(oI,4);
% 
% %%Get Centers
% oProps = regionprops(lI,oI,'MinorAxisLength','MajorAxisLength','Area',...
%     'WeightedCentroid');
% Cent = I * 0;
% for o = 1:oNum
%    Cent(round(oProps(o).WeightedCentroid(2)),round(oProps(o).WeightedCentroid(1))) = 1; 
% end
% 
% Cent = imdilate(Cent,strel('disk',3));
% colI = uint8(fitH(oI));
% colI(:,:,2) = fitH(Cent);
% colI(:,:,3) = fitH(tpI);
% 
% image(colI)
% 
% %% Get Edges within Wells
% % [lI oNum] = bwlabel(oI>1,8); %%label objects
% % eI = I*0;
% % image(fitH(eI>0))
% % for o = 1: oNum
% % o
% %     [oy ox] = find(lI == o);
% %     sI = dI(min(oy):max(oy),min(ox):max(ox));
% %     sI(sI>0) = sI(sI>0) - min(sI(sI>0));
% %     min(oy)
% %     idO = sub2ind(size(sI),oy-min(oy)+1,ox-min(ox)+1);
% %     esI = edge(sI,'Canny',[.001 .1],3);
% %     lE = bwlabel(esI,8);
% %     goodEdges = lE(idO);
% %     goodEdges = unique(goodEdges(goodEdges>0));
% %     badEdges = setdiff(1:max(lE(:)),goodEdges);
% %     subplot(1,2,1)
% %     image(fitH(sI))
% %     subplot(1,2,2)
% %     image(fitH(lE))
% %     pause
% %     for e = 1:length(badEdges)
% %         lE(lE==badEdges(e))=0;
% %     end
% %     eI(min(oy):max(oy),min(ox):max(ox))= (lE>0) * o;
% % 
% % end
% % image(fitH(eI)),colormap colorcube(255),pause(.1)
% % %% combine objects with edges
% % 
% % eIa = eI * 0;
% % eProps = regionprops(eI,'Orientation','Centroid','ConvexHull');
% % for i = 1:length(eProps)
% %     cPoly = eProps(i).ConvexHull;
% %     eIa(poly2mask(cPoly(:,1),cPoly(:,2),ys,xs)) = 1;
% % end
% % image(eIa*1000),pause(.1)
% % colI(:,:,3) = fitH(eIa);
% % colI(:,:,2) = fitH(I);
% % colI(:,:,1) = fitH(tpI>0);
% % image(colI),pause(.1)
% % 
% % [lI oNum] = bwlabel(eIa>0,4); %%label objects
% %% select sections
% % 
% % rProps = regionprops(lI,cI2,'Area','EulerNumber','Orientation',...
% %     'WeightedCentroid','MeanIntensity','MajorAxisLength',...
% %     'MinorAxisLength','Perimeter','Extent','Extrema','ConvexHull',...
% %     'PixelIdxList','BoundingBox');
% % 
% % minAx = [rProps.MinorAxisLength];
% % majAx = [rProps.MajorAxisLength];
% % 
% % [IDX,C,sumd,D] = kmeans(minAx,2);
% % 
% % scatter(minAx(IDX == 1),majAx(IDX == 1),'.','r');
% % hold on
% % scatter(minAx(IDX == 2), majAx(IDX == 2),'.','b');
% % hold off
% % 
% % %%Remove small
% % useK = find(IDX == find(C == min(C)));
% % for i = 1: length(useK)
% %     lI(lI == useK(i)) = 0;
% % end
% % 
% % colI(:,:,3) = fitH(lI>0);
% % image(colI),pause(.1)
% % 
% % [lI oNum] = bwlabel(lI>0,4);
% % rProps = regionprops(lI,'Centroid','ConvexHull');
% % 
% % centI = I*0;
% % for i = 1:length(rProps)
% %     centI(round(rProps(i).Centroid(2)),round(rProps(i).Centroid(1))) = 1;
% % end
% % 
% % centI = imdilate(centI,strel('disk',3));
% % colI(:,:,2) = centI*1000;
% % 
% % image(colI),pause(.01)
% % 
% 
% 
% toc