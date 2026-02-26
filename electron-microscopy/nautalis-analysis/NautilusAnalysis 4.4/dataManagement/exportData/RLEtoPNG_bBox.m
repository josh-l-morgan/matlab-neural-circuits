% function[rleDat] = getRLEtoSubs(TPN,miplevel);
%
% %Please connect to VAST with vasttools first!
%
% global vdata;
%
% tic

if ~exist('vdata','var')
    vasttools
end
global vdata

if ~vdata.state.isconnected
    vdata.vast.connect('127.0.0.1',22081,1000)
end

TPN = GetMyDir;
TPN2 = [TPN(1:end-1) '_lowRes\'];
mkdir(TPN2)
writeLowRes = 1;

%%

zStep = 16; %100 for mip3, 20 for mip 1
%xyStep = 4096;
xStep = 16*16*16;
yStep = 16*16*16;
miplevel = 0;
checkMip = 2;
checkDown = 2^checkMip;
%% Get seg data
disp('Getting all segment data for targeted fetching')
nr = vdata.vast.getnumberofsegments;
bBox = zeros(nr-1,6);
for i = 1:nr-1
    data = vdata.vast.getsegmentdata(i);
    bBox(i,:) = double(data.boundingbox);
end
max(bBox,[],1)
useBox = sum(bBox>=0,2)==6;
bBox = bBox(useBox,:);
allBBox = [min(bBox(:,1:3),[],1)'  max(bBox(:,4:6),[],1)']+1;
allBBox2 = allBBox/checkDown
allBBox2(:,1) = fix(allBBox2(:,1));
allBBox2(:,2) = ceil(allBBox2(:,2));



blockDim = [yStep xStep zStep];



% (36673, 26090, 687)

info=vdata.vast.getinfo();
segNum = vdata.vast.getnumberofsegments
segData = vdata.vast.getallsegmentdata;
%segNames = vdata.vast.getallsegmentnames


minx=0; miny=0; minz=0;
maxx=double(bitshift(info.datasizex,-miplevel));
maxy=double(bitshift(info.datasizey,-miplevel));
%same as: maxy=floor(info.datasizey/2^miplevel);
maxz=double(info.datasizez);
vSize = [maxy-miny+1 maxx-minx+1 maxz-minz+1];
surfonlyflag=0;


datSize = [maxx; maxy; maxz];
zSteps = 1:zStep:maxz;
xSteps = 1:xStep:maxx;
ySteps = 1:yStep:maxy;

% zSteps = 687:zStep:687+200;
% xSteps = 36673:xStep:36673 + 10000;
% ySteps = 26090:yStep:26090 + 10000;

startTime = clock

i = 0;
numBlocks = length(zSteps) * length(ySteps) ;
countYvox = 0;
countYseg = 0;
startYtime = tic;

for z = 1:length(zSteps);
    disp(sprintf('fetching sections %d-%d of %d',zSteps(z),zSteps(z)+zStep-1,maxz))
    for y = 1:length(ySteps);
        yDur = toc(startYtime);
        disp(sprintf('%0.2f Mseg at %.2f Mvox/sec, %.2f Mseg/sec',countYseg,countYvox/yDur,...
            countYseg/yDur))
        startYtime = tic;
        countYvox = 0;
        countYseg = 0;

        for x = 1:length(xSteps);
            
            if ~mod(y,100)
                disp(sprintf('running %d of %d',y,length(ySteps)))
            end
            
            
            %%Get region
            getReg = [xSteps(x)-1 min(xSteps(x)+xStep -2,allBBox(1,2));...
                ySteps(y)-1 min(ySteps(y)+yStep -2,allBBox(2,2));
                zSteps(z)-1 min(zSteps(z)+zStep -2,allBBox(3,2))];
            getBBox = [getReg(:,1)' getReg(:,2)'];
            hit = bBoxIntersect(getBBox,bBox);
            if hit
            regSize = getReg(:,2)-getReg(:,1) + 1;
            countYvox = countYvox + prod(regSize)/1000000;
            
            %%Test for data
            getReg2 = getReg;
            datSize2 = allBBox;
            datSize2(1:2) = datSize2(1:2)/(2^checkMip);
            datSize2 = round(datSize);
            getReg2(1:2,1:2) = round(getReg2(1:2,1:2)/2^checkMip);
            getReg2(1:2,1) = getReg2(1:2,1) - 2;
            getReg2(1:2,2) = getReg2(1:2,2) + 2;
            getReg2(getReg2(:,1)<0,1) = 0;
            getReg2(getReg2(:,2)>datSize2,2) = datSize2(getReg2(:,2)>datSize2);
            regSize2 = getReg2(:,2)-getReg2(:,1) + 1;
            Iv2 =  zeros(regSize2(1),regSize2(2),regSize2(3));
            [rleimg2, res]=vdata.vast.getsegimageRLE(checkMip,...
                getReg2(1,1),getReg2(1,2), getReg2(2,1),getReg2(2,2),...
                getReg2(3,1),getReg2(3,2),surfonlyflag);
            
            if writeLowRes
                lastX = 0;
                oIDs = double(rleimg2(1:2:end));
                rls = double(rleimg2(2:2:end));
                for r = 1:length(rls)
                    runID = lastX+1:lastX+rls(r);
                    
                    if oIDs(r)>0
                        Iv2(runID) = oIDs(r);
                        countYseg = countYseg+length(runID)/1000000;
                    end
                    lastX = lastX + rls(r);
                end
                
                %%write 2D images
                iSums = squeeze(sum(sum(Iv2)));
                writeIs = find(iSums);
                for zi = 1:length(writeIs)
                    writeI = writeIs(zi);
                    imName = sprintf('tileLowRes_r%d_c%d_s%d.png',y,x,...
                        getReg(3,1)-1 + writeI);
                    imwrite(Iv2(:,:,writeI),[TPN2 imName],'BitDepth',16)
                end
            end
            %toc
            
            %%Get full resolution
            if sum(iSums) %if low res found something
                Iv = zeros(regSize(1),regSize(2),regSize(3),'uint16');
                [rleimg, res]=vdata.vast.getsegimageRLE(miplevel,...
                    getReg(1,1),getReg(1,2), getReg(2,1),getReg(2,2),...
                    getReg(3,1),getReg(3,2),surfonlyflag);
                
                %%write run length into 3D image
                lastX = 0;
                oIDs = double(rleimg(1:2:end));
                rls = double(rleimg(2:2:end));
                for r = 1:length(rls)
                    runID = lastX+1:lastX+rls(r);
                    
                    if oIDs(r)>0
                        Iv(runID) = oIDs(r);
                    end
                    lastX = lastX + rls(r);
                end
                
                %%write 2D images
                iSums = squeeze(sum(sum(Iv)));
                writeIs = find(iSums);
                for zi = 1:length(writeIs)
                    writeI = writeIs(zi);
                    imName = sprintf('tile_r%d_c%d_s%d.png',y,x,...
                        getReg(3,1)-1 + writeI);
                    imwrite(Iv(:,:,writeI),[TPN imName],'BitDepth',16)
                end
                
            end %if low res found something
            else % if in seg range
                %disp('no segmentation within fetch range')
            end
        end
    end
end









 %             %%Get whole image
                %             tic
                %             [img, res]=vdata.vast.getsegimageraw(miplevel,...
                %                 getReg(1,1),getReg(1,2), getReg(2,1),getReg(2,2),...
                %                 getReg(3,1),getReg(3,2));
                %             Iv(:) = img;
                %             toc
