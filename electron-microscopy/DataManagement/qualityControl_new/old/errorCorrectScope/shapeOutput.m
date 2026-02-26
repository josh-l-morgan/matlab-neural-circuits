%function[retake] = shapeOutput(wif,focSec)
if ~exist('wif','var')
    wif= GetMyWafer;
end

if ~exist('focSec','var')
    load([wif.dir 'focSec.mat'])
end

if exist([wif.dir 'userEval.mat'])
    load([wif.dir 'userEval.mat'])
else
    userEval = 0
end

if isfield(userEval,'thresh')
    thresh = userEval.thresh;
    qualityThresh = thresh(1);
else
    qualityThresh = 20;
end

saturationThresh = 5;
evalThresh = 3; %minimum user Evaluation

retake.qualityThresh = qualityThresh;
retake.saturationThresh = saturationThresh;
retake.evalThresh = evalThresh;
retake.time = datestr(clock);

rNum = 0;
numSec = length(wif.sec);
for s = 1:numSec

    [tree, rootname, dom]=xml_read(wif.sec(s).xml);


    numTiles = length(wif.sec(s).tile);
    for t = 1:numTiles
        quality = focSec(s).tile(t).quality;
        percentSat = focSec(s).tile(t).percentSaturation;
        if isfield(userEval,'qualVals')
            qualSec = userEval.qualVals(:,:,s);
            rating = qualSec(t);
            bad = rating< evalThresh;
        else
            bad = quality< qualityThresh;
            bad = bad + (percentSat > saturationThresh);
            rating = 'none';
        end
        if isfield(focSec(s).tile(t),'retake')
            bad = focSec(s).tile(t).retake;
        end
%% if bad
        if bad
            rNum = rNum + 1;

            if ~exist('tree','var')
                [tree, rootname, dom]=xml_read(wif.sec(s).xml);
            end

            %% Find focus target
            TargetX = tree.Tiles.Tile(t).TargetStageX;
            TargetY = tree.Tiles.Tile(t).TargetStageY;
            FOV = tree.TileInfo.FOV.CONTENT;
            ScanRot = tree.ReferenceInfo.Beam.ScanRot;
            focTarg = focSec(s).tile(t).focusTarget;

            shiftY = -FOV * (focTarg(1)-.5);
            shiftX = -FOV * (focTarg(2)-.5);
            [r h ] = cart2pol(shiftX,shiftY);
            newR = r + ScanRot/360*2*pi;
            [newX newY] = pol2cart(newR,h);
            
            focTargYum = TargetY +  newY ;
            focTargXum = TargetX + newX ;            
            
            
%pause
            %% Record

            tile(rNum).imTargX = TargetX / (-1 * 10^6);
            tile(rNum).imTargY = TargetY / (-1 * 10^6);
            tile(rNum).focTargX = focTargXum / (-1 * 10^6);
            tile(rNum).focTargY = focTargYum / (-1 * 10^6);
            tile(rNum).startWD = tree.ReferenceInfo.Beam.WD;
            tile(rNum).startSX = tree.ReferenceInfo.Beam.StigX;
            tile(rNum).startSY = tree.ReferenceInfo.Beam.StigY;
            tile(rNum).ScanRotation = ScanRot;
            tile(rNum).StageRotation = tree.ReferenceInfo.Stage.Rot;
            tile(rNum).TileWidth = tree.TileInfo.TileWidth;
            tile(rNum).TileHeight = tree.TileInfo.TileHeight;
            tile(rNum).FOV = FOV;
            tile(rNum).DwellTime = tree.DwellTime.CONTENT/1000;
            tile(rNum).Brightness = tree.Tiles.Tile(t).Brightness;
            tile(rNum).Contrast = tree.Tiles.Tile(t).Contrast;
            tile(rNum).PixelSize = tree.PixelSize.CONTENT * 10^-9;
            
            
            tileInfo(rNum).path = wif.sec(s).tile{t};
            nam = wif.sec(s).tile{t};
            slash = find(nam == '\');
            tileInfo(rNum).name = nam(slash(end)+1:end);
            tileInfo(rNum).sectionNum = s;
            tileInfo(rNum).rc = tree.Tiles.Tile(t).ATTRIBUTE;
            
            oldQuality(rNum).quality = quality;
            oldQuality(rNum).percentSaturation = percentSat;
            oldQuality(rNum).userEval = rating;
            
            %%Find focus mag
            
            
        end %if bad


    end

end %end tiles
clear tree

if exist('tile','var')
retake.tiles = tile;
retake.tileInfo = tileInfo;
retake.oldQuality = oldQuality;

else
    retake = 0;

end
save([wif.dir 'retake.mat'],'retake')



