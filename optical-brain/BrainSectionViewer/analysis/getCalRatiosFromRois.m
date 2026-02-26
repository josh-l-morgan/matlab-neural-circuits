

SPN = 'E:\SeanMcCracken\testCaNoise\395 (processed)\'
clear vReport

dSPN = dir([SPN '*.tif'])

iNams = {dSPN.name};        %name of file needs to be idential to its ROI name for characters after roiPreTag 

showAll = 1;

roiPreTag = 'RoiSet_';

for i = 1:length(iNams)
    
    %% Read image
    clear I
    nam = iNams{i};
    
    info = imfinfo([SPN nam]);
    L = length(info);
    ys = info(1).Height;
    xs = info(1).Width;
    zs = L/2;
    
    Iraw = zeros(ys,xs,zs);
    for p = 1:L
        Iraw(:,:,p) = imread([SPN nam],p);
    end
    
    zeroMode = (2^16)/2;%mode(Iraw(:));  % correct for unsigned read
   
    %%make color image
    I{1} = Iraw(:,:,1:2:end-1);
    I{2} = Iraw(:,:,2:2:end);
    
    %% Filter data
    if 0
        
        %% Filter data
        colormap gray(255)
        for c = 1:2
            Ic = I{c};
            I1 = imgaussfilt3(Ic,[2 2 3]); %center filter
            %I2 = imgaussfilt3(Ic,[7 7 12]); %surround filter
            %I3 = I1-I2;
            If{c} = I1;
        end
        I = If;
  
    end
    
     %%normalize to background
    zeroMode1 = mode(I{1});
    zeroMode2 = mode(I{2});  
     
    a = I{1}== zeroMode1;
    for b = 1:size(a,3);
        b
        image(a(:,:,b)*1000)
        drawnow
    end
    
    I{1} = I{1}-zeroMode1;
    I{2} = I{2}-zeroMode2;
    
    
    %% show color
    Icol = zeros(ys,xs,3,zs,'double');
    for  p = 1:zs
        Icol(:,:,1,p) = I{1}(:,:,p);
        Icol(:,:,2,p) = I{2}(:,:,p);
        Icol(:,:,3,p) = I{1}(:,:,p);
    end
    Icol = Icol - min(Icol(:));
    Icol = Icol * 255/max(Icol(:));
    Icol = uint8(Icol);
    for p = 1:zs
        image(Icol(:,:,:,p));
        pause(.1)
    end
    
    %% Read rois
    
    roiNam = [roiPreTag nam(1:end-4) '.zip'];
    roiFileName = [SPN roiNam];
    
    [sROI] = ReadImageJROI(roiFileName)
    rYs = sROI{1}.vnRectBounds(2);
    rXs = sROI{1}.vnRectBounds(4);
    [sRegions] = ROIs2Regions(sROI, [ys xs]);
    idxs = sRegions.PixelIdxList;
    rs = length(idxs);
    
    mask = zeros(ys,xs,zs);
    clear val
    bin = 1;
    for r = 1 :length(idxs)
        sroi = sROI{r};
        roiZ = sroi.vnPosition(2);
        
        [yt xt] = ind2sub([ys xs],idxs{r});
        
        if ~isempty(yt)
        
        newInd1 = sub2ind([ys xs zs],xt,yt,repmat(roiZ,[length(yt) 1]));
        newInd2 = sub2ind([ys xs zs],xt,yt,repmat(roiZ,[length(yt+1) 1]));
        newInd3 = sub2ind([ys xs zs],xt,yt,repmat(roiZ,[length(yt-1) 1]));
        newInd4 = sub2ind([ys xs zs],xt,yt,repmat(roiZ,[length(yt+2) 1]));
        newInd5 = sub2ind([ys xs zs],xt,yt,repmat(roiZ,[length(yt-12) 1]));
        
        
        newInd = [newInd1; newInd2; newInd3; newInd4; newInd5];

        
        [newY newX newZ] = ind2sub([ys xs zs],newInd);
        
        val{1}{r} = I{1}(newInd);
        val{2}{r} = I{2}(newInd);
        
        if showAll
            mask = mask*0;
            mask(newInd) = 1;
            showI = Icol(:,:,:,roiZ) * 4;
            showI(:,:,3) = uint8(mask(:,:,roiZ)*100);
            image(showI)
            pause(.01)
            drawnow
        end
        else
            val{1}{r} = [];
            val{2}{r} = [];
        end
    end
    
    clear vCount vStd rS rM vErr vRatio vMean 
    zVal = 1.645;% 1.645 = 95%, 1.282 corresponds to 90%
    
    for r = 1:rs
        
        v1 = val{1}{r};
        v2 = val{2}{r};
        
        vCount(r,1) = length(v1);
        vCount(r,2) = length(v2);
        vStd(r,1) = std(v1);
        vStd(r,2) = std(v2);
        vMean(r,1) = mean(v1);
        vMean(r,2) = mean(v2);
        vRatio(r,1) = mean(v1)/mean(v2);
        
    end
    
    vSE = vStd./sqrt(vCount);
    vHigh = vMean + vSE * zVal;
    vLow = vMean - vSE * zVal;
    vErr(:,1) = vHigh(:,1) ./ vLow(:,2);
    vErr(:,2) = vLow(:,1) ./ vHigh(:,2);
    
    [sortR idx] = sort(vRatio,'ascend');
    plot(vRatio(idx),'k');
    hold on
    plot(vErr(idx,:),'r')
    hold off
    pause(1)
    
    vReport.header = {'Ratio', 'High error','Low error','Mean 1','Mean 2','Pix count','STD1', 'STD2'};
    vReport.imageName{i} = nam;
    vReport.dat{i} = cat(2,vRatio,vErr,vMean,vCount,vStd);
    vReport.zeroModes = [zeroMode1 zeroMode2];
    
end

save([SPN 'vReport.mat'],'vReport');















