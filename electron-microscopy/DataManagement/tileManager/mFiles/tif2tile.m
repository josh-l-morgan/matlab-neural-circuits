function[tileName] = tif2tile(tifName)
        
tilePos = regexp(tifName,'Tile_');
tifName = tifName(tilePos(1):end);

unds = regexp(tifName,'_');
if length(unds) == 3
    tpos = regexp(tifName,'.tif');
    tileName = tifName(1:tpos(1)-1);
elseif length(unds) == 4
            tileName = tifName(1:unds(4)-1);
else
    tpos = regexp(tifName,'.tif');
    tileName = tifName(1:tpos(1)-1);
    disp(sprintf('Wrong number of underscores.'))
end
