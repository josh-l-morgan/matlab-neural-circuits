function[allVox] = bridgeConMat(allVox)


clf
subs = dropSubs(allVox.subs);
showSteps = 0;

seed = allVox.seedInd;

renderConK(subs,[],[1 0 0])

allVox.conMatOld = allVox.conMat;
bSubs = dropSubs(allVox.subs);
rangeA = max(bSubs,[],1);
allInds = sub2ind(rangeA,bSubs(:,1),bSubs(:,2),bSubs(:,3));

con = allVox.conMat;
deep = sum(con>0,2);


segCount = 0;
pairs = [];
clf





%% Segment by connectivity
segID = zeros(size(con,1),1);
lastSegID = segID;
segCount = 0
clf
hold on
while 1
    segCount = segCount + 1;
    disp(sprintf('filling segment %d',segCount));
    seed = find(segID==0,1);
    segID(seed) = segCount;
    
    %%Spread
    while 1
        for d = 1:26
            moveID =  (segID>0) & (con(:,d)>0);
            segID(con(moveID,d)) = segID(moveID);
        end
        %                 image(segID*100)
        %                 drawnow
        if sum(abs(lastSegID-segID))==0
            break
        end
        lastSegID =segID;
    end
    if sum(segID==segCount) == 0
        break
    else
        renderConK(subs(segID==segCount,:),[],[rand rand rand]),drawnow,
    end
end


uSegs = unique(segID);
segC = hist(segID,uSegs);
smallSeg = find(segC == min(segC),1);
if showSteps
    clf
    for i = 1:length(uSegs);
        renderConK(subs(segID==uSegs(i),:),[],[rand rand rand]),drawnow
    end
end


while 1
    
    uSegs = unique(segID);
    disp(sprintf('current number of segments = %d',length(uSegs)))
    if length(uSegs)<2
        break
    end
    segC = hist(segID,uSegs);
    smallSeg = uSegs(find(segC == min(segC),1));
    
    %%Find nearest
    isFilled = ~(segID == smallSeg);
    isNot = segID==smallSeg;
    fullPosF = subs(isFilled,:);
    fullPosN = subs(isNot,:);
    clf
    if showSteps, renderConK(fullPosF,[],[1 1 1],.1), end
  
    if ~sum(isNot)
        break
    end
    if showSteps, renderConK(fullPosN,[],[0 0 1],.1),pause(.01), end
    
    %% Find potential overlaps
    
    boxSize = 10;
    rangeF = [min(fullPosF,[],1)-10;max(fullPosF,[],1)+10];
    rangeN = [min(fullPosN,[],1)-10; max(fullPosN,[],1)+10];
    rangeBuf = boxSize/2; %extension of boundaries to include voxels as close, increases with itterations
    rangeBuf2 = 0;
    buffObj = 10; %buffer around object boundaries

    %% DSpos
    dsPos = 5;
    fullPosFds = fullPosF;
    fullPosFds = round(fullPosFds /dsPos)+1;
    rangeFds = max(fullPosFds,[],1);
    fullPosFdsInd = sub2ind(rangeFds,fullPosFds(:,1),fullPosFds(:,2),fullPosFds(:,3));
    fullPosFdsInd = unique(fullPosFdsInd);
    [y x z] = ind2sub(rangeFds,fullPosFdsInd);
    fullPosFds = ([y x z] - 1) * dsPos;
    
    fullPosNds = fullPosN;
    fullPosNds = round(fullPosNds /dsPos)+1;
    rangeNds = max(fullPosNds,[],1);
    fullPosNdsInd = sub2ind(rangeNds,fullPosNds(:,1),fullPosNds(:,2),fullPosNds(:,3));
    fullPosNdsInd = unique(fullPosNdsInd);
    [y x z] = ind2sub(rangeNds,fullPosNdsInd);
    fullPosNds = ([y x z] - 1) * dsPos;    
    
    
    for b = 1:20
        
        closeN = [];
        closeF = [];
        %% Do boxes overlap?
        for y = floor(rangeN(1,1)/boxSize):ceil(rangeN(2,1)/boxSize)
            
            disp(sprintf('checking %d of %d',y,ceil(rangeN(2,1)/boxSize)))
            
            rangeBuf2 = rangeBuf; %!!!!!!!!!!
            
            for x = floor(range(1,2)/boxSize):ceil(rangeN(2,2)/boxSize)

                for z = floor(range(1,3)/boxSize):ceil(rangeN(2,3)/boxSize)
                    
                    
                     hasF =  ((fullPosFds(:,1)+rangeBuf)>((y-1)*boxSize)) & ((fullPosFds(:,1)-rangeBuf)<(y*boxSize)) & ...
                        ((fullPosFds(:,2)+rangeBuf)>((x-1)*boxSize)) & ((fullPosFds(:,2)-rangeBuf)<(x*boxSize)) & ...
                        ((fullPosFds(:,3)+rangeBuf)>((z-1)*boxSize)) & ((fullPosFds(:,3)-rangeBuf)<(z*boxSize));
                    
                     hasN =  ((fullPosNds(:,1)+rangeBuf2)>((y-1)*boxSize)) & ((fullPosNds(:,1)-rangeBuf2)<(y*boxSize)) & ...
                        ((fullPosNds(:,2)+rangeBuf2)>((x-1)*boxSize)) & ((fullPosNds(:,2)-rangeBuf2)<(x*boxSize)) & ...
                        ((fullPosNds(:,3)+rangeBuf2)>((z-1)*boxSize)) & ((fullPosNds(:,3)-rangeBuf2)<(z*boxSize));
                    hasBoth = sum(hasF) & sum(hasN);
                    
                    
                    if hasBoth % if both filled and unfilled downsampled voxels are present in range
                    
                    
                    inRange =  ((subs(:,1)+rangeBuf)>((y-1)*boxSize)) & ((subs(:,1)-rangeBuf)<(y*boxSize)) & ...
                        ((subs(:,2)+rangeBuf)>((x-1)*boxSize)) & ((subs(:,2)-rangeBuf)<(x*boxSize)) & ...
                        ((subs(:,3)+rangeBuf)>((z-1)*boxSize)) & ((subs(:,3)-rangeBuf)<(z*boxSize));
                    
                    
                    posF = subs(isFilled & inRange,:);
                    posN = subs(isNot & inRange,:);
                    
                    if ~isempty(posF) & ~isempty(posN) %Search range relative to contents of box
                        
                        
                        if 0
                        subRangeF = [min(posF,[],1)-buffObj;max(posF,[],1)+buffObj];
                        subRangeN = [min(posN,[],1)-buffObj; max(posN,[],1)+buffObj];
                        
                        
                        
                        inRangeF =  ((posF(:,1)+rangeBuf)>((subRangeN(1,1)))) & ((posF(:,1)-rangeBuf)<(subRangeN(2,1))) & ...
                            ((posF(:,2)+rangeBuf)>((subRangeN(1,2)))) & ((posF(:,2)-rangeBuf)<(subRangeN(2,2))) & ...
                            ((posF(:,3)+rangeBuf)>((subRangeN(1,3)-1))) & ((posF(:,3)-rangeBuf)<(subRangeN(2,3)));
                        
                        inRangeN =  ((posN(:,1)+rangeBuf)>((subRangeF(1,1)))) & ((posN(:,1)-rangeBuf)<(subRangeF(2,1))) & ...
                            ((posN(:,2)+rangeBuf)>((subRangeF(1,2)))) & ((posN(:,2)-rangeBuf)<(subRangeF(2,2))) & ...
                            ((posN(:,3)+rangeBuf)>((subRangeF(1,3)-1))) & ((posN(:,3)-rangeBuf)<(subRangeF(2,3)));
                        
                        closeF = cat(1,closeF,posF(inRangeF,:));
                        closeN = cat(1,closeN,posN(inRangeN,:));
                        else
                            closeF = cat(1,closeF,posF);
                            closeN = cat(1,closeN,posN);
                        end
                        
                            if showSteps
                                'bark'
                                clf
                                renderConK(fullPosN,[],[0 0 1],.1),pause(.01), 
                                hold on
                                if ~isempty(posF)
                                renderConK(posF,[],[1 0 0],.2),pause(.01), 
                                end
                                if ~isempty(posN)
                                  renderConK(posN,[],[0 1 0],.2),pause(.01), 
                                end
                               drawnow
                            end

                        
                    end %is pos
                    
                    end % is both in down sampled
                end %run z
            end %run x
        end %runy
        
        if ~isempty(closeF) & ~isempty(closeN)
            
            uF = unique(sub2ind(rangeA,closeF(:,1),closeF(:,2),closeF(:,3)));
            uN = unique(sub2ind(rangeA,closeN(:,1),closeN(:,2),closeN(:,3)));
            [y1 x1 z1] = ind2sub(rangeA,uF);
            closeF = [y1 x1 z1];
            [y2 x2 z2] = ind2sub(rangeA,uN);
            closeN = [y2 x2 z2];
            
            if showSteps, renderCon(closeF+.5,[],[1 0 0]), end
            if showSteps, renderCon(closeN+.5,[],[0 1 0]), end
            pause(.01)
            break
            
        end
        
        rangeBuf = rangeBuf * 2;
        
    end
    
    
    %% 
       
    %%%%%%%%dists = sqrt((posF(:,1)-posN(:,1)').^2 + (posF(:,2)-posN(:,2)').^2 + (posF(:,3)-posN(:,3)').^2);
    %dists = sqrt((closeF(:,1)-closeN(:,1)').^2 + (closeF(:,2)-closeN(:,2)').^2 + (closeF(:,3)-closeN(:,3)').^2);
    
    
    
    dists = zeros(size(closeF,1),size(closeN,1),'single');
    for d = 1:size(closeN,1)
        dists(:,d) = sqrt((closeF(:,1) - closeN(d,1)).^2 + (closeF(:,2) - closeN(d,2)).^2 + ...
            (closeF(:,3) - closeN(d,3)).^2);
    end
    
%     
%     tic
%     minDist = inf;
%     for d = 1:size(closeN,1)
%         distsD = sqrt((closeF(:,1) - closeN(d,1)).^2 + (closeF(:,2) - closeN(d,2)).^2 + ...
%             (closeF(:,3) - closeN(d,3)).^2);
%         minD = min(distsD);
%         if minD<minDist
%             n = find(distsD == minD,1);
%             f = d;
%         end
%     end
%     toc
    
    
    
    
    minDist = min(dists(:));
    [f n] = find(dists == minDist,1);
    
    hitF = closeF(f,:);
    hitN = closeN(n,:);
    hitFind = sub2ind(rangeA,hitF(1),hitF(2),hitF(3));
    hitNind = sub2ind(rangeA,hitN(1),hitN(2),hitN(3));
       
    pF = find((allInds==hitFind) & isFilled,1);
    if showSteps
        renderCon(subs(pF,:)+.5,[],[1 1 0],1);
    end
    pN =  find((allInds==hitNind) & ~isFilled);
    if showSteps
        renderCon(uniqueSubs(subs(pN,:))+.5,[],[0 1 1],1);
    end
    pause(.01);
    pN = setdiff(pN,pF); % in case of multiple hits
    pN = pN(1);
    pairs  = cat(1,pairs,[pF pN]);
    
    seed = pN; %ready for next spreading
    
    %%change small seg to connected found seg ID
    foundSeg = segID(pF);
    if showSteps
        clf
        renderConK(subs,[],[0 0 1]),
        renderConK(subs(segID==smallSeg,:),[],[0 1 0]),
        renderConK(subs(segID==foundSeg,:),[],[1 0 0]),
        drawnow
    end
    segID(segID==smallSeg) = foundSeg;

    
    if showSteps
        hold on
        for b = 1:size(pairs,1)
            nList = pairs(b,:);
            runE = bSubs(nList,:);
            plot3(runE(:,1),runE(:,2),runE(:,3),'linewidth',3,'color','w')
        end
        
        pause(.01)
    end
    
end

%% change conMat

for p = 1:size(pairs,1)
    
    con1 = con(pairs(p,1),:);
    con2 = con(pairs(p,2),:);
    
    empty1 = find(con1==0,1);
    empty2 = find(con2==0,1);
    
    con(pairs(p,1),empty1) = pairs(p,2);
    con(pairs(p,2),empty2) = pairs(p,1);
    
end

allVox.conMat = con;
allVox.voxBridges = pairs; 

if 1 %% show
    clf
    renderConK(subs,allVox.conMat,[1 0 0])
    hold on
        for b = 1:size(pairs,1)
            nList = pairs(b,:);
            runE = bSubs(nList,:);
            plot3(runE(:,1),runE(:,2),runE(:,3),'linewidth',3,'color','w')
        end  
    
    ax = gca;
    ax.Clipping = 'off';
    
end

pause(.01)








