
function[obV] = rle2subs(rleimg,param,segNum)


%% convert to voxel list

bufVnum = 100000;
clear obV 
obV{segNum} = [];
trackO = zeros(segNum,1);
    
vSize = [        param.maxx - param.minx + 1 ...
        param.maxy - param.miny + 1 ...
        param.maxz - param.minz + 1];
    
    lastPos = 0;
    
    for r = 2:2:length(rleimg)
        oID = rleimg(r-1);
        if oID>0
            linIdx = lastPos + 1: lastPos + double(rleimg(r));
            [x y z] = ind2sub( vSize, linIdx);
            y = y + double(param.miny) - 1;
            x = x + double(param.minx) - 1;
            z = z + double(param.minz) - 1;
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
            lastPos = lastPos + double(rleimg(r));
        end
        
    end
    
    
    for t = 1:segNum %ditch buffer
        obV{t} = obV{t}(1:trackO(t),:);
    end

% 
% targ = find(trackO == max(trackO),1);
% subs = obV{targ}(:,1:2);
% subs(:,1) = subs(:,1)-min(subs(:,1))+1;
% subs(:,2) = subs(:,2) - min(subs(:,2))+1;
% maxI = max(subs,[],1);
% showI = zeros(maxI);
% showI(sub2ind(maxI,subs(:,1),subs(:,2))) = 1000;
% image(showI)
% 
% 





