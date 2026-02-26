SPN = 'D:\LGNs1\rawMontages\';
WaferList = {'w030','w031','w032','w033','w034','w035','w036'}

TPN = 'D:\LGNs1\isStack3-2\';

if ~exist(TPN,'dir');
    mkdir(TPN);
end


dSPN = dir(SPN); dSPN = dSPN(3:end);
dNams = cat(2,{dSPN.name});
isDirectory = cat(1,dSPN.isdir);
dNams = dNams(isDirectory);

getTileName = 'Tile_r3-c2';

montageDir = regexp(dNams,'Montage');
montageDir = ~cellfun('isempty',montageDir);

sortNams = {};
for w = 1:length(WaferList)
   
    waferName = WaferList{w};
    rightWafer = regexp(dNams,waferName);
    rightWafer = ~cellfun('isempty',rightWafer);
    rightSections =dNams(rightWafer& montageDir);
    sortSections = sort(rightSections);
    wafSecs{w} = sortSections;
     
end


c = 0;
for w = 1:length(wafSecs);
    secNam = wafSecs{w};
    for s = 1:length(secNam);
        disp(sprintf('Downsampling wafer %d of %d section %d of %d',...
            w,length(wafSecs),s,length(secNam)))
        nam = secNam{s};
        secDir  = dir([SPN nam]); secDir = secDir(3:end);
        tileNams = cat(1,{secDir.name});
        tileNams = sort(tileNams);
        findTile = regexp(tileNams,getTileName);
        findTile = ~cellfun('isempty',findTile);
        findTif = regexp(tileNams,'.tif');
        findTif = ~cellfun('isempty',findTif);
        targTif = find(findTile&findTif,1,'first');
        foundName = tileNams{targTif};
        if exist([SPN nam '\' foundName],'file')
            
            
            
            %%
            
            FileTif=[SPN nam '\' foundName];
InfoImage=imfinfo(FileTif);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
FinalImage=zeros(nImage,mImage,NumberImages,'uint16');
FileID = tifflib('open',FileTif,'r');
rps = tifflib('getField',FileID,Tiff.TagID.RowsPerStrip);
 
for i=1:NumberImages
   tifflib('setDirectory',FileID,i);
   % Go through each strip of data.
   rps = min(rps,mImage);
   for r = 1:rps:mImage
      row_inds = r:min(mImage,r+rps-1);
      stripNum = tifflib('computeStrip',FileID,r);
      FinalImage(row_inds,:,i) = tifflib('readEncodedStrip',FileID,stripNum);
   end
end
tifflib('close',FileID);
            
            
            
            
            
            
            
            
            
            %%
            
        tic
            I = imread([SPN nam '\' foundName]);
        toc
        tic
        dsI = imresize(I,1/8,'nearest');
        toc
        pause
        c = c+1;
         newName = [TPN 'ds' zerobuf(c,5) '_' foundName];
         imwrite(dsI,newName,'Compression','none')
        end
    end
end