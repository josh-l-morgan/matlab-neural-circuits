
TPN = GetMyDir;

load([TPN 'rleDat.mat'])



%% convert to voxel list

bufVnum = 100000;
segNum = rleDat.segNum;
blockNum = length(rleDat.rleimg);
clear obV stackO
%obV{segNum} = [];



for i = 1:blockNum
    trackO = zeros(segNum,1);
    obV{segNum} = [];
    
    tic
    disp(sprintf('running block %d of %d',i,blockNum))
    rleimg = double(rleDat.rleimg{i});
    
    vSize = [        rleDat.param(i).maxx - rleDat.param(i).minx + 1 ...
        rleDat.param(i).maxy - rleDat.param(i).miny + 1 ...
        rleDat.param(i).maxz - rleDat.param(i).minz + 1];
    
    lastPos = 0;
    
    for r = 2:2:length(rleimg)
        oID = rleimg(r-1);
        if oID>0
            linIdx = lastPos + 1: lastPos + rleimg(r);
            [x y z] = ind2sub( vSize, linIdx);
            y = y + rleDat.param(i).miny - 1;
            x = x + rleDat.param(i).minx - 1;
            z = z + double(rleDat.param(i).minz) - 1;
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
    toc
    
    
    for t = 1:segNum %ditch buffer
        obV{t} = obV{t}(1:trackO(t),:);
    end
    
    stackO{i}.obV = obV;
end



targ = find(trackO == max(trackO),1);
subs = obV{targ}(:,1:2);
subs(:,1) = subs(:,1)-min(subs(:,1))+1;
subs(:,2) = subs(:,2) - min(subs(:,2))+1;
maxI = max(subs,[],1);
showI = zeros(maxI);
showI(sub2ind(maxI,subs(:,1),subs(:,2))) = 1000;
image(showI)







