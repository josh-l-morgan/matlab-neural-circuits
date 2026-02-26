function[] = downSampWithReplacement(SPN,TPN,secNam,s,getTileName,c);

nam = secNam{s};
secDir  = dir([SPN nam]); secDir = secDir(3:end);
tileNams = cat(1,{secDir.name});
tileNams = sort(tileNams);
findTile = regexp(tileNams,getTileName);
findTile = ~cellfun('isempty',findTile);
findTif = regexp(tileNams,'.tif');
findTif = ~cellfun('isempty',findTif);
targTif = find(findTile&findTif,1,'first');
if ~isempty(targTif)
    foundName = tileNams{targTif};
    newName = [TPN 'ds' zeroBuf(c,5) '_' foundName];
    
    if exist([SPN nam '\' foundName],'file') & ~exist(newName)
        
        
        I = imread([SPN nam '\' foundName]);
        
       shift1 = zeros(size(I,1)*size(I,2)/8,8);
        
        
        
        
        imwrite(dsI,newName,'Compression','none')
    end
end


