
chan(1).dir = 'D:\LGNs1\HP_processing\rayAlign01\Vast\export01_notMem\'
chan(2).dir =  'D:\LGNs1\HP_processing\rayAlign01\Vast\export01_mem\'
chan(3).dir = 'D:\LGNs1\HP_processing\rayAlign01\images\'

TPN = 'D:\LGNs1\HP_processing\rayAlign01\Vast\export01_col\';
if ~exist(TPN,'dir')
    mkdir(TPN)
end


for c = 1:length(cDir)
    iList = {};
    dirC = dir(chan(c).dir);
    for i = 1:length(dirC);
        nam = dirC(i).name;
        if sum(regexp(nam,'.tif')) | sum(regexp(nam,'.png'))
           iList{length(iList)+1} = nam; 
        end
    end
    chan(c).iList = iList; 
end


scaling = [ 1000 1000 1];
for i = 1:length(chan(1).iList)
   for c = 1:length(chan)
      colI(:,:,c) = uint8(imread([chan(c).dir chan(c).iList{i}])* scaling(c));
   end   
   iName = sprintf('memYN_%05.0f.tif',i);
   imwrite(colI,[TPN iName]);
end
