
SPN = GetMyDir;
TPN = [SPN(1:end-1) '_ds\'];
mkdir(TPN)

%% read names
dSPN = dir(SPN);
iNam = {};
for i  =   1:length(dSPN)
   nam = dSPN(i).name;
   if sum(regexp(nam,'.tif'))
      iNam{length(iNam)+1} = nam; 
   end    
end


%%

for i = 1:length(iNam);
    disp(sprintf('%d of %d',i,length(iNam)));
    I = imread([SPN iNam{i}]);
    
    dsI = imresize(I,.25,'bicubic');
    
    imwrite(dsI,[TPN iNam{i}]);
end