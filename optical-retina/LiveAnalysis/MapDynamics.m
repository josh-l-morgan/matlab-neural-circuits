

%    regNames = {}
%    comGo = []

for i = 1: 6%size(regNames,1)   %% Run all cells
   TX = load([regNames{i,1} regNames{i,2}]);
   Points{i}=[mean([TX(:,17),TX(:,23)],2) mean([TX(:,18),TX(:,24)],2)];
    
   movers = comGo(comGo(:,1)== i,2:4);
   
   %% Draw 
   Point = Points{i};
   Point = round(Point);
   Point(Point<1)=1;
   imSize = max(Point);
   pointInd = sub2ind(imSize,Point(:,1),Point(:,2));
   allP = zeros(imSize); newP = allP; goneP = allP;
   allP(pointInd) = 200;
   
   formed = find(movers(:,2));
   lost = find(movers(:,3));
   
   newP(pointInd(movers(formed,1))) = 200;
   goneP(pointInd(movers(lost,1))) = 200;
   
   I = cat(3,goneP,newP,allP);
   SE = strel('disk',5);
   for c = 1: 3
      I(:,:,c) = imdilate(I(:,:,c), SE ); 
   end
   
   images{i} = uint8(I);
end


c = 1
while 1
   c = c + 1; 
   image(images{mod(c,6)+1}) 
   pause
end

