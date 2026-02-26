function[regDat]  = registerBrainSections(regDat)

%%Give the function regDat containint I1 (fixed image), I2 (Image to move)
%%ifShow (1 shows images), dsamp (amount to downsamp to speed up, 2 is
%%default, ms1 and ms2 for mexican hat filter, maxIterations default = 100

if isfield(regDat,'ifShow')
    ifShow = regDat.ifShow;
else
    ifShow = 1;
end

if ~isfield(regDat, 'dsamp')
    regDat.dsamp = 2;
end
dsamp = regDat.dsamp;

if isfield(regDat,'maxIterations')
    maxIterations = regDat.maxIterations;
else
    maxIterations = 100;
end

if ~isfield(regDat, 'mergeChannels')
    regDat.mergeChannels = 0;
end
mergeChannels = regDat.mergeChannels;


if isfield(regDat,'app')
        ax = regDat.app.ax;
elseif ifShow
    ax = gca;
end


ifCutEdges = 1;

if isfield(regDat,'ms1')
    ifFilt = 1;
    ms1 = regDat.ms1;
    ms2 = regDat.ms2;
    g1 = fspecial('gaussian',ms2*4,ms1);
    g1 = g1/sum(g1(:));
    g2 = fspecial('gaussian',ms2*4,ms2);
    g2 = g2/sum(g2(:));
    kern = g1-g2;
else 
    ifFilt = 0;
end



I1 = double(regDat.I1);
I2 = double(regDat.I2);

if mergeChannels
    I1 = mean(I1,3);
    I2 = mean(I2,3);
end


mI1 = mean(I1,3);
mI2 = mean(I2,3);

if ifShow
    title(ax,'before alignment')
    imagePair(mI2, mI1,ax);
    pause(.1)
end

numChan = size(I1,3);

[optimizer, metric] = imregconfig('multimodal');
 optimizer.InitialRadius = 0.0001;
% optimizer.Epsilon = 1.5e-8;
% optimizer.GrowthFactor = 1.001;
optimizer.MaximumIterations = maxIterations;

clear tform
T = zeros(3,3,numChan);
cc = zeros(1,1,numChan);
trans = zeros(numChan,2);
rot = zeros(numChan,1);
cc = zeros(numChan,1);
for c = 1: size(I1,3)
    I1c = I1(:,:,c);
    I2c = I2(:,:,c);
    I1cDs = imresize(I1c,1/dsamp);
    I2cDs = imresize(I2c,1/dsamp);
    if ifFilt
        I1cDs = imfilter(I1cDs,kern,'same');
        I2cDs = imfilter(I2cDs,kern,'same');
        cut = ceil(size(kern,1)/2);
        I1cDs = I1cDs(cut+1:end-cut,cut+1:end-cut);
        I2cDs = I2cDs(cut+1:end-cut,cut+1:end-cut);
    end

    tform = imregtform(I2cDs ,I1cDs ,'rigid',optimizer,metric);
    T(:,:,c) = tform.T;
    I2r = imwarp(I2cDs ,tform,'OutputView',imref2d(size(I1cDs)));
    
    rot(c) = tform.RotationAngle;
    trans(c,:) = tform.Translation;
    coefs = corrcoef(I1cDs,I2r);
    cc(c) = coefs(1,2);

    if ifShow
        imagePair(I2r, I1cDs,ax);
        title(ax,sprintf('channel %d, cc = %0.2f',c,cc(c)))
        drawnow
    end
end

mTrans = median(trans,1);
mRot = median(rot);
mRot = deg2rad(mRot);
mT = [cos(mRot) sin(mRot) 0;...
    -sin(mRot) cos(mRot) 0;...
    mTrans(1) mTrans(2) 1];

%mT = median(T,3);
%mT(1:2,1:2) = [1 0;0 1];
%mt(3,3) = 1;
%mT = sum(T .* cc,3)/sum(cc)
mT(3,1:2) = mT(3,1:2) * dsamp;
tform.T = mT;
mI2r = imwarp(mI2,tform,'OutputView',imref2d(size(mI1)));

if ifShow
    title(ax,'with alignment')
    imagePair(mI2r, mI1,ax);
end
regDat.tform = tform;









