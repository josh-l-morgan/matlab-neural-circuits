function[tforms] = alignMatSectionsFunc(regDat)
%% Align images.  Images must be matfiles in folder
SPN = regDat.SPN;
global globBB
% regDat.dsamp = 2;
% regDat.mergeChannels = 1; %else channels will be aligned separately and median value used
% regDat.rLim = 30; %maximum rotation is 30 degrees
% regDat.shiftLim = .3; % maximum shift is 30% (0.3) of image size
% regDat.ms1 = 1;
% regDat.ms2 = 4;
ax = regDat.app.ax;
flipSec = globBB.show.flip;

groupDir = SPN;
if exist(groupDir,'dir')
    gDir = dir([groupDir 'sec*.mat']);
    nams = {gDir.name};
    startI = round(length(nams)/2);
    tforms{startI} = affine2d;

    %%forward step alignment
    load([groupDir nams{startI}],'I');
    if flipSec(startI)
        I = fliplr(I);
    end
    regDat.I1 = I;

    mIlast = mean(I,3);
    for z = startI+1:length(nams);
        textStr = (sprintf('running forward align %d of %d, group %d',z,length(nams),i));
        regDat.app.textOut.Value = textStr;
        mI1 = mean(regDat.I1,3);
        load([groupDir nams{z}],'I');
        if flipSec(z)
            I = fliplr(I);
        end
        regDat.I2 = I;
        mI2 = mean(regDat.I2,3);


        regDat = registerBrainSections(regDat);

        %clear transform if out of bounds
        shiftLim = regDat.shiftLim .* size(mI1,1);
        r  = abs(acosd(regDat.tform.T(1,1)));
        maxShift = max(abs(regDat.tform.T(3,[1 2])));
        if ~((r <= regDat.rLim) & (maxShift <= shiftLim))
            regDat.tform.T = [1 0 0; 0 1 0 ; 0 0 1];
        end

        tforms{z} = regDat.tform;


        %%transform
        It = imwarp(regDat.I2,tforms{z},'OutputView',imref2d(size(mI1)));
        mIt = mean(It,3);
        title(ax,'after alignment')
        imagePair(mI1, mIt,ax);
        pause(.1)
        %Swich current to previous
        regDat.I1 = It;
    end

    %%back step alignment
    load([groupDir nams{startI}],'I');
    if flipSec(startI)
            I = fliplr(I);
    end
    regDat.I1 = I;
    mIlast = mean(I,3);
    for z = startI-1:-1:1;
        textStr = (sprintf('running backwards align %d to 1, group %d',z,i));
        regDat.app.textOut.Value = textStr;
        mI1 = mean(regDat.I1,3);
        mI2 = mean(regDat.I2,3);
        load([groupDir nams{z}],'I');
        if flipSec(z)
            I = fliplr(I);
        end
        regDat.I2 = I;
        regDat = registerBrainSections(regDat);

        %clear transform if out of bounds
        shiftLim = regDat.shiftLim .* size(mI1,1);
        r  = abs(acosd(regDat.tform.T(1,1)));
        maxShift = max(abs(regDat.tform.T(3,[1 2])));
        if ~((r <= regDat.rLim) & (maxShift <= shiftLim))
            regDat.tform.T = [1 0 0; 0 1 0 ; 0 0 1];
        end

        tforms{z} = regDat.tform;

        %%transform
        It = imwarp(regDat.I2,tforms{z},'OutputView',imref2d(size(mI1)));
        mIt = mean(It,3);
        title(ax,'after alignment')
        imagePair(mI1, mIt,ax);
        pause(.1)

        %Swich current to previous
        regDat.I1 = It;
    end

    save([groupDir 'tforms.mat'],'tforms');
end

