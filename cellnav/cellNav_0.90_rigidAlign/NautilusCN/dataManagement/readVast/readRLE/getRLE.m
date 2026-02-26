function[rleDat] = getRLE(TPN,miplevel);

%Please connect to VAST with vasttools first!

global vdata;

tic

if ~exist('miplevel','var')
    miplevel=3;
end
zStep = 64; %100 for mip3, 20 for mip 1
xyStep = 1024;



%TPN = GetMyDir;
info=vdata.vast.getinfo();
segNum = vdata.vast.getnumberofsegments
segData = vdata.vast.getallsegmentdata;
%segNames = vdata.vast.getallsegmentnames

rleDat.info = info;
rleDat.segNum = segNum;
rleDat.segData = segData;
%rleDat.segNames = segNames;

minx=0; miny=0; minz=0;
maxx=bitshift(info.datasizex,-miplevel);
maxy=bitshift(info.datasizey,-miplevel);
%same as: maxy=floor(info.datasizey/2^miplevel);
maxz=info.datasizez;
vSize = [maxy-miny+1 maxx-minx+1 maxz-minz+1];
surfonlyflag=0;



zSteps = 1:zStep:maxz; 
%zSteps = floor(maxz/2); %TEMP for debug
xSteps = 1:xyStep:maxx;
ySteps = 1:xyStep:maxy;

startTime = clock

i = 0;
for z = 1:length(zSteps);
    startZ = zSteps(z);
    stopZ = min(startZ+zStep-1,maxz);
                disp(sprintf('running %d of %d blocks',z,length(zSteps)))

    for y = 1:length(ySteps);
        startY = ySteps(y);
        stopY = min(startY+xyStep-1,maxy);
        
        for x = 1:length(xSteps);
            
            i = i + 1;
            
            startX = xSteps(x);
            stopX = min(startX+xyStep-1,maxx);
            
            
            [rleimg, res]=vdata.vast.getsegimageRLE(miplevel,startX-1,stopX-1,startY-1,stopY-1,startZ-1,stopZ-1,surfonlyflag);
            rleDat.rleimg{i} = rleimg;
            
            
            param.stepSize(i) = stopZ-startZ+1;
             param.miplevel = miplevel;
%             param.minx = minx;
%             param.maxx = maxx;
%             param.miny = miny;
%             param.maxy = maxy;
            param.minx = startX;
            param.maxx = stopX ;
            param.miny = startY;
            param.maxy = stopY;
            param.minz = startZ;
            param.maxz = stopZ;
            param.surfonlyflag = surfonlyflag;
            rleDat.param(i) = param;
        end
        
    end
end
        save([TPN 'rleDat.mat'],'rleDat','-v7.3')
        
        
        stopTime = clock
