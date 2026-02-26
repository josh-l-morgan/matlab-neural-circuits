function[retake] = shapeOutput(wif,focSec)
if ~exist('wif','var')
    wif= GetMyWafer;
end
if ~exist('focSec','var')
    load([wif.dir 'focSec.mat'])

end

qualityThresh = 20;
saturationThresh = 5;
evalThresh = 3; %minimum user Evaluation

retake.qualityThresh = qualityThresh;
retake.saturationThresh = saturationThresh;
retake.evalThresh = evalThresh;
retake.time = datestr(clock);


rNum = 0;
numSec = length(wif.sec)
for s = 1:numSec

    [tree, rootname, dom]=xml_read(wif.sec(s).xml);


    numTiles = length(wif.sec(s).tile);
    for t = 1:numTiles
        quality = focSec(s).tile(t).quality;
        percentSat = focSec(s).tile(t).percentSaturation;
        if isfield(focSec(s).tile(t),'userEval')
            userEval = focSec(s).tile(t).userEval;
            bad = userEval< evalThresh;
        else
            bad = quality< qualityThresh;
            bad = bad + (percentSat > saturationThresh);
            userEval = 'none';
        end

        if bad
            rNum = rNum + 1;

            if ~exist('tree','var')
                [tree, rootname, dom]=xml_read(wif.sec(s).xml);
            end

            %% Find focus target
            TargetX = tree.Tiles.Tile(t).TargetStageX;
            TargetY = tree.Tiles.Tile(t).TargetStageY;
            FOV = tree.TileInfo.FOV.CONTENT;
            focTarg = focSec(s).tile(t).focusTarget;
            focTargYum = TargetY - FOV/2 + FOV * focTarg(1);
            focTargXum = TargetX - FOV/2 + FOV * focTarg(2);

            %% Record

            tile(rNum).imTargX = TargetX;
            tile(rNum).imTargY = TargetY;
            tile(rNum).focTargX = focTargXum;
            tile(rNum).focTargY = focTargYum;
            tile(rNum).startWD = tree.Tiles.Tile(t).WD;
            tile(rNum).startSX = 0;
            tile(rNum).startSY = 0;
            tile(rNum).TileWidth = tree.TileInfo.TileWidth;
            tile(rNum).TileHeight = tree.TileInfo.TileHeight;
            tile(rNum).FOV = FOV;
            tile(rNum).DwellTime = tree.DwellTime.CONTENT;

            tileInfo(rNum).name = wif.sec(s).tile{t};
            tileInfo(rNum).sectionNum = s;
            tileInfo(rNum).rc = tree.Tiles.Tile(t).ATTRIBUTE;
            
            oldQuality(rNum).quality = quality;
            oldQuality(rNum).percentSaturation = percentSat;
            oldQuality(rNum).userEval = userEval;
        end


    end

end %end tiles
clear tree


retake.tiles = tile;
retake.tileInfo = tileInfo;
retake.oldQuality = oldQuality;

save([wif.dir 'retake.mat'],'retake')





