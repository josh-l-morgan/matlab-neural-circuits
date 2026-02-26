%Please connect to VAST with VastTools first!

global vdata;

TPN = GetMyDir;
info=vdata.vast.getinfo();
%%? seg = vdata.vast.GETSEGMENTDATA()
%%? vdata.vast.GETNUMBEROFSEGMENTS

nr=vdata.vast.getnumberofsegments()

miplevel=3;
minx=0; miny=0; minz=900;
maxx=bitshift(info.datasizex,-miplevel)-1;
maxy=bitshift(info.datasizey,-miplevel)-1;
%same as: maxy=floor(info.datasizey/2^miplevel);
maxz=info.datasizez;
vSize = [maxy-miny+1 maxx-minx+1 maxz-minz+1];
surfonlyflag=0;

zStep = 100;
zSteps = 1:zStep:maxz;

startTime = clock
for i = 1:length(zSteps);
    disp(sprintf('running %d of %d blocks',i,length(zSteps)))
    tic
    startZ = zSteps(i);
    if i ==  length(zSteps)
        stopZ = maxz;
    else
        stopZ = zSteps(i+1)-1;
    end
    [rleimg, res]=vdata.vast.getsegimageRLE(miplevel,minx,maxx,miny,maxy,startZ,stopZ,surfonlyflag);
  
    allRLEi{i} = rleimg;
    stepSize(i) = stopZ-startZ+1;
    toc
end

rleDat.stepSize = stepSize;
rleDat.allRLEi = allRLEi;

save([TPN 'rleDat.mat'],'rleDat')


stopTime = clock


%% convert to voxel list

bufObNum = nr;

bufVnum = 10000;
trackO = zeros(bufObNum,1); 
clear obV
obV{bufObNum} = [];

for i = 1: length(allRLEi)
    disp(sprintf('running block %d of %d',i,length(allRLEi)))
    rleimg = double(allRLEi{i});
    lastPos = 0;
    for r = 2:2:length(rleimg)
        oID = rleimg(r-1);
        if oID>0
            linIdx = lastPos + 1: lastPos + rleimg(r);
            [y x z] = ind2sub( vSize, linIdx);   
            lastPos = linIdx(end);
            
            %%Find current and projected voxel list length
            endV = trackO(oID) + length(y);
            [currN threes] = size(obV{oID});
            
            if endV>currN % Add more voxel slots ( to avoid constant growth)
                obV{oID}(endV+bufVnum,1:3) = [0 0 0];
            end
            
            obV{oID}(trackO(oID)+1:endV,:) = [y' x' z'];
            
            trackO(oID) = endV; %update number of real voxels for object
        else
            lastPos = lastPos + rleimg(r); 
        end
        
    end    
    
end









