

%%
Xum = 6000 * 1000
Yum = 3000 * 1000
Zum = 7000 * 1000

cutComp = .8;
secThick = 40
fillVolume = 0.5236

pixSize = 10
6
secOrder = [ 2 3 1]; % section dim, compression dim, width

volDim = [Yum Xum Zum];
volDim = volDim(secOrder);

sectionNumber = volDim(1)/ secThick
cuttingLength = sectionNumber * volDim(2) * fillVolume / 1000000000

voxelNumber = volDim(1)/secThick * volDim(2)/pixSize * volDim(3)/pixSize * fillVolume





