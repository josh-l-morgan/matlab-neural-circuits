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
        replaceWith = mean(I(:));
        %replaceWith = mean([replaceWith tooHigh]);
       shift1 = zeros(8,size(I,1)*size(I,2)/8);
       shift1(:) = I(:); 
       tooHigh = 230;
       highI = shift1>=tooHigh;
       shift1(highI) = replaceWith;
       %numHigh = sum(highI,1);
       sumAll = sum(double(shift1),1);
       testShift = zeros(3200,25600);
       testShift(:) = sumAll(:);
      % image(testShift/8);
       
       shift2 = zeros(3200,8,length(sumAll)/3200/8);
       %shift2a = shift2;
       shift2(:) = sumAll(:);
       %shift2a(:) = numHigh(:);
       %finalHigh = squeeze(sum(shift2a,2));
       finalSum = squeeze(sum(shift2,2));
       
       dsMeansSc = finalSum./64;
       %dsMeansSc(finalHigh>=63) = replaceWith;
       image(256-dsMeansSc)
        
        imwrite(dsI,newName,'Compression','none')
    end
end


