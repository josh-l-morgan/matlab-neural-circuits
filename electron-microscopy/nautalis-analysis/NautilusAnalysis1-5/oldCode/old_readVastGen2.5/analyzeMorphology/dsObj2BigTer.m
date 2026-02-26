
%% load data
clear all
MPN = 'D:\LGNs1\Segmentation\VAST\S8\joshm\export+14+04+27_mat\'
TPN = [MPN 'skel\'];
load([MPN 'dsObj.mat'])
load([MPN 'obI.mat'])

minTouch = 100;

%% load paths
mainDir = pwd;
slash = regexp(mainDir,'\');
mainDir = mainDir(1:slash(end));
path(path,sprintf('%sanalyzeMorphology\\skeletonizeByShortestDistance',mainDir))
path(path,sprintf('%sanalyzeMorphology',mainDir))
path(path,sprintf('%sviewCells',mainDir))
path(path,sprintf('%sanalyzeColors',mainDir))
path(path,sprintf('%sanalyzeConnectivity',mainDir))
path(path,sprintf('%svast2ob',mainDir))




%% Define variables
lookDist = 10;
allSubs = cat(1,dsObj.subs);
maxSub = max(allSubs,[],1)+10;
clear allSubs


%% downSample object by 8

dSamp = 4;

%% down sample

for i = 1:length(dsObj)
    
    disp(sprintf('%d of %d',i,length(dsObj)));
    sub = ceil(double(dsObj(i).subs)/dSamp);
    sub(sub<1) = 1; %%!!! Problem getting subs should never get zeros
    maxSub = max(sub,[],1);
    if ~isempty(sub)
        inds = sub2ind(maxSub,sub(:,1),sub(:,2),sub(:,3));
        uInds = unique(inds);
        if length(uInds>1)
            hInds = hist(inds,uInds);
        else
            hInds = length(inds)
            'bark'
        end
        [y x z] = ind2sub(maxSub,uInds);
        dsDsObj(i).subs = cat(2,uint16(y),uint16(x),uint16(z));
        %dsObj(i).n = uint16(hInds);
    end
    
end

%save([MPN 'dsObj.mat'],'dsObj','-v7.3')
clear vastSubs




%% make ball filter
lookDist = 12;
ballSize = [lookDist*2+1  lookDist*2+1 lookDist*2+1];
ball = strel('ball',lookDist,lookDist);
ball = ones(ballSize);
ballInd = find(ball);
[bY bX bZ] = ind2sub(ballSize,ballInd);
dists = sqrt((bY - lookDist - 1).^2 + (bX - lookDist - 1).^2 + (bZ - lookDist - 1).^2);
ball(dists>lookDist) = 0;
%ball(dists<lookDist-1.8) = 0;
colormap gray(256)

image(sum(ball,3)*10)

%%
for post = 1: length(postList)
    post
    mid2 = getCellSubs(obI,dsDsObj,postList(post));
    mid2(:,1:2) = mid2(:,1:2) + 1024/8;
    %% lesmid
    maxMid = max(mid2,[],1);
    midInd = sub2ind(maxMid,mid2(:,1),mid2(:,2),mid2(:,3));
    midInd = unique(midInd);
    mid3 = zeros(length(midInd),3);
    [mid3(:,1) mid3(:,2) mid3(:,3)] = ind2sub(maxMid,midInd);
    conMat = obj2con(mid3);
    getSurf = sum(conMat>0,2)<23;
    mid4 = mid3(getSurf,:);
    
    tic
    diSub2 = dilateSub(mid4,ball);
    toc
    
    
    diSub2(diSub2(:,1)>maxSub(1),1) = maxSub(1);
    diSub2(diSub2(:,2)>maxSub(1),2) = maxSub(2);
    diSub2(diSub2(:,3)>maxSub(1),3) = maxSub(3);
    
    allPostDiSub{post} = diSub2;
    
end




