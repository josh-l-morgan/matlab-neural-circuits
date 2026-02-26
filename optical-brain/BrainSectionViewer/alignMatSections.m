%% Align images.  Images must be matfiles in folder
SPN = SPN;
clf
clear regDat
regDat.dsamp = 2;
regDat.mergeChannels = 1; %else channels will be aligned separately and median value used
regDat.rLim = 30; %maximum rotation is 30 degrees
regDat.shiftLim = .3; % maximum shift is 30% (0.3) of image size
% regDat.ms1 = 1;
% regDat.ms2 = 4;
af = figure;
for i = 1:length(imageTags)
    groupDir = [SPN imageTags{i} '\'];
    if exist(groupDir,'dir')
        if length(nams)>1
            gDir = dir([groupDir 'sec*.mat']);
            nams = {gDir.name};
            startI = round(length(nams)/2);
            tforms{startI} = affine2d;

            %%forward step alignment
            load([groupDir nams{startI}],'I');
            regDat.I1 = I;
            mIlast = mean(I,3);
            for z = startI+1:length(nams);
                disp(sprintf('running forward align %d of %d, group %d',z,length(nams),i))
                mI1 = mean(regDat.I1,3);
                load([groupDir nams{z}],'I');
                regDat.I2 = I;
                mI2 = mean(regDat.I2,3);
                regDat = registerBrainSections(regDat);
               
                %clear transform if out of bounds
                r  = abs(acosd(regDat.tform.T(1,1)));
                maxShift = max(abs(regDat.tForm.T(3,[1 2])));
                if ~((r <= regDat.rLim) & (maxShift <= regDat.shiftLim))
                    regDat.tform.T = [1 0 0; 0 1 0 ; 0 0 1];
                end
                
                tforms{z} = regDat.tform;


                %%transform
                It = imwarp(regDat.I2,tforms{z},'OutputView',imref2d(size(mI1)));
                mIt = mean(It,3);
                imshowpair(mI1, mIt,'Scaling','joint')
                title('after alignment')
                pause(1)
                %Swich current to previous
                regDat.I1 = It;
            end

            %%back step alignment
            load([groupDir nams{startI}],'I');
            regDat.I1 = I;
            mIlast = mean(I,3);
            for z = startI-1:-1:1;
                disp(sprintf('running backwards align %d to 1, group %d',z,i))
                mI1 = mean(regDat.I1,3);
                mI2 = mean(regDat.I2,3);
                load([groupDir nams{z}],'I');
                regDat.I2 = I;
                regDat = registerBrainSections(regDat);
                tforms{z} = regDat.tform;

                %%transform
                It = imwarp(regDat.I2,tforms{z},'OutputView',imref2d(size(mI1)));
                mIt = mean(It,3);
                imagePair(mI1, mIt);
                title('after alignment')
                pause(1)

                %Swich current to previous
                regDat.I1 = It;
            end

            save([groupDir 'tforms.mat'],'tforms');
        end
    end
end

pause(2)
close(af)