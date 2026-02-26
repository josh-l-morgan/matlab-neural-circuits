
SPN = GetMyDir
TPN = GetMyDir


%%

zStart = 1;
zStop = 200;
xStart = 6488;
xStop = 9280;
yStart = 6203;
yStop = 7732;


iNams = {};
dSPN = dir(SPN);
for i = 1:length(dSPN)
   nam = dSPN(i).name;
   if sum(regexp(nam,'.tif'));
       iNams{length(iNams)+1} = nam;
   end
    
end
%%

for i = 1:min(zStop,length(iNams));
    i
   I = imread([SPN iNams{i}],'PixelRegion',{[yStart yStop] [xStart xStop]});
    imwrite(I,[TPN iNams{i}]);
    
end