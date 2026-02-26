function[]=axDotNeuriteFinder(TPN,DPN) %#ok<INUSL>
% 
% TPN = GetMyDir;
% DPN=[TPN 'I\'] %use when running as stand alone script
%% Dot Finder 
%F3 has been modified to handle 8bit single tiff stacks
%large files will be broken into smaller blocks and then recombined

%% Get file names
tic

%Get directory name
f=find(DPN=='\');
f2=f(size(f,2)-1);
f3=f(size(f,2)-2);
TPN=DPN(1:f2); %Define target folder (one level up from files)
if isdir('./history')==0, mkdir('./history'); end %create directory to store steps
save(['./history' TPN(f3:f2-1)],'TPN') %record path in history folder

colormap gray(255) %standard grey colormap
if isdir([TPN 'temp'])==0, mkdir([TPN 'temp']); end %create directory to store steps
if isdir([TPN 'data'])==0, mkdir([TPN 'data']); end %create directory to store steps
if isdir([TPN 'pics'])==0, mkdir([TPN 'pics']); end %create directory to store steps
 


%% Enter Variables
%%Image Variables
% xyum=.103;
% zum=.3;
% % aspect=zum/xyum;% ratio of z to xy dimentions
% channels=3;
BlockSize=500; %aproximate size of individual processing blocks
BlockBuffer=30; %amount of block to be cut off as edge buffer

%%Dot criteria
step=2; %Sensitivity=grey value step of iterative threshold (2 was standard)
MaxDot=7^3;  %Maximum Dot Volume (in pixels)= maximum dot size for iterative threshold (6^3 was standard)
MinDot=6;  %Minimum Dot Volume (in pixels)= minimum dot size for iterative threshold   (10 was standard)
% PunctaThreshold=1; %Minimum Number of Steps Passed = centroids less than puncta threshold are zeroed. 
% EdgeOfPeak=.5; %Determine Dot Edge = ratio of edge brightness to peak brightness (.5 was standard)
% minFilledVolume=3; %minimum number of pixels in final object contour (20 was standard)
% RoundThreshold=50; %minimum roundness threshold (60 was standard)


%% READ IMAGE
%%Find Images
d=dir(DPN); d=d(3:size(d,1)); %get number of files in directory
Rzs = size(d,1); %find number of planes
I(:,:,:)=imread([DPN d(round(Rzs/2)).name]);
[Rys,Rxs,dummy]=size(I);   %#ok<NASGU> %Read first pic 
IgRaw=zeros(Rys,Rxs,Rzs,'uint8');

%%Figure out channels
if size(I,3)==1, 
    puncta=1; %if only one channel
elseif size(I,3)==2
    puncta=1; %if only two channels
    neurite=2;
elseif sum(sum(I(:,:,3)))==0
    puncta=1; %if third channel is blank
    neurite=2;
else
    puncta=2; %if third channel is not blank
    neurite=1;
end

for i=1:Rzs
    I(:,:,:)=imread([DPN d(i).name]);
    iPunctaRaw(:,:,i)=uint8(I(:,:,puncta)); %ColorSeperate
    iNeuriteRaw(:,:,i)=uint8(I(:,:,neurite)); %ColorSeperate
end
clear I
save([TPN 'temp\iRaw'],'iPunctaRaw', 'iNeuriteRaw')

%% Find Image Stats
if ~exist('iPunctaRaw','var') %load green
    load([TPN 'temp\iRaw']);
end
iPunctaSignal = iPunctaRaw(:);
iPunctaSignal(iPunctaSignal == 0) = [];
IgStat = sort(iPunctaSignal);
Gmode = IgStat(round(0.75*length(IgStat)));
clear IgStat

%% SubSample
if Rxs>BlockSize
    NumBx=round(Rxs/BlockSize);
    Bxc=fix(Rxs/NumBx); 
else
    Bxc=Rxs;
    NumBx=1;
end

if Rys>BlockSize
    NumBy=round(Rys/BlockSize);
    Byc=fix(Rys/NumBy); 
else
    Byc=Rys; 
    NumBy=1; 
end

if Rzs>BlockSize
    NumBz=round(Rzs/BlockSize);
    Bzc=fix(Rzs/NumBz); 
else
    Bzc=Rzs; 
    NumBz=1; 
end

%% Run Blocks
for Bz=1:NumBz, for By=1:NumBy, for Bx=1:NumBx      %#ok<ALIGN>

PercentBlocksDone = ((Bz-1)*NumBy*NumBx+Bz   +  ...
    (By-1) * NumBx + By  + Bx)/(NumBz*NumBx*NumBy)  %#ok<NASGU,NOPRT>

%Find real territory
Tystart=(By-1)*Byc+1;
Txstart=(Bx-1)*Bxc+1;
Tzstart=(Bz-1)*Bzc+1;
if By<Byc, Tyend=By*Byc; else Tyend=Rys; end
if Bx<Bxc, Txend=Bx*Bxc; else Txend=Rxs; end
if Bz<Bzc, Tzend=Bz*Bzc; else Tzend=Rzs; end

%Find buffered Borders (extend to image boarders for last blocks in row and column)
yStart=Tystart-BlockBuffer;
yStart(yStart<1)=1;
yEnd=Tyend+BlockBuffer;
yEnd(yEnd>Rys)=Rys;
xStart=Txstart-BlockBuffer;
xStart(xStart<1)=1;
xEnd=Txend+BlockBuffer;
xEnd(xEnd>Rxs)=Rxs;
zStart=Tzstart-BlockBuffer;
zStart(zStart<1)=1;
zEnd=Tzend+BlockBuffer;
zEnd(zEnd>Rzs)=Rzs;


%%Subsample Green Channel
if ~exist('iPunctaRaw','var')%load green if necessary
    load([TPN 'temp\iRaw'])
end 
iPuncta = single(iPunctaRaw(yStart:yEnd,xStart:xEnd,zStart:zEnd)); 
iNeurite = single(iNeuriteRaw(yStart:yEnd,xStart:xEnd,zStart:zEnd));
clear iPunctaRaw iNeuriteRaw
[ys,xs,zs]=size(iPuncta);

%% Median filter 
iPunctaMed=zeros(ys,xs,zs,'single');
iNeuriteMed=zeros(ys,xs,zs,'single');
for i=1:zs        
    iPunctaMed(:,:,i)=medfilt2(iPuncta(:,:,i),[3,3]);
    iNeuriteMed(:,:,i)=medfilt2(iNeurite(:,:,i),[3,3]);
end
clear iPuncta iNeurite


%% FIND DOTS Green Channel%%
peakMap = zeros(ys,xs,zs,'uint8');  %set up matrix to map peaks
thresholdMap = zeros(ys,xs,zs,'uint8');   %set up matrix to sum passed thresholds
% indVect = 1:ys*xs*zs;
maxIntensity = uint8(max(iPunctaMed(:)));

for i = maxIntensity:-step:Gmode
    %run thresholds through all relevant intensities
    clear Igl labels
    [Igl,labels] = bwlabeln(iPunctaMed>i,6);%label each area to check

    %reduce bitdepth if possible
    if labels<65536
        Igl=uint16(Igl);
    end
    if labels <= 1
        labels = 2;
    end
    nPixel = hist(Igl(Igl>0), 1:labels);
    %run all lables
    for p=1:labels
        % Morphology Filter, Puncta size criteria
        pixelIndex = find(Igl==p);
        if nPixel(p) < MaxDot && nPixel(p) > MinDot
            %identify peak in labeled object
            if sum(peakMap(pixelIndex))== 0
                peakValue = max(iPunctaMed(pixelIndex));
                peakIndex = find(Igl==p & iPunctaMed==peakValue);
                if numel(peakIndex) > 1
                    peakIndex = peakIndex(round(numel(peakIndex)/2));
                end
                [y,x,z] = ind2sub([ys xs zs], peakIndex);
                %Register in peak map
                peakMap(y,x,z) = 1;
            end
        else
            Igl(pixelIndex)=0;
        end
    end
    %%Add all passing labeled objects to thresholdMap
    thresholdMap(Igl>0)=thresholdMap(Igl>0)+1;
end 
disp('iterative threshold done')
clear Igl peakIndex


%% FIND DOT CONTOUR AND DIVIDE IF MULTIPLE PEAKS WITHIN
[thresholdLabel nLabels] = bwlabeln(thresholdMap, 6);

if nLabels<65536 %reduce bit depth if possible
    thresholdLabel = uint16(thresholdLabel);
end

for i = 1:nLabels
    peakIndex = find(thresholdLabel==i & peakMap>0);
    thresholdPeak = thresholdMap(peakIndex);
    nPeaks = numel(peakIndex);
    [yPeak xPeak zPeak] = ind2sub([ys xs zs], peakIndex);
    if nPeaks == 1
        if  (yPeak-1 < BlockBuffer/2 && yStart > 1) ||...
                (xPeak-1 < BlockBuffer/2 && xStart > 1) ||...
                (zPeak-1 < BlockBuffer/2 && zStart > 1) ||...
                (ys-yPeak < BlockBuffer/2 && yEnd < Rys) ||...
                (xs-xPeak < BlockBuffer/2 && xEnd < Rxs) ||...
                (zs-zPeak < BlockBuffer/2 && zEnd < Rzs)
        else
            cutOff = 0.5 * thresholdPeak;
            contourIndex = find(thresholdLabel==i & thresholdMap>=cutOff);
            [yContour xContour zContour] = ind2sub([ys xs zs],contourIndex);
            if ~exist('Dots','var')
                lastEntry = 0;
            else
                lastEntry = size(Dots.Pos,1);  
            end
            next = lastEntry+1;
            Dots.Pos(next,:) = [yPeak+yStart-1, xPeak+xStart-1,...
                zPeak+zStart-1];
            Dots.Vox(next).Pos = [yContour+yStart-1,...
                xContour+xStart-1, zContour+zStart-1];
            Dots.Vox(next).Ind = sub2ind([Rys Rxs Rzs],...
                Dots.Vox(next).Pos(:,1), Dots.Vox(next).Pos(:,2),...
                Dots.Vox(next).Pos(:,3));
            Dots.Vol(next) = numel(contourIndex);
            Dots.ITMax(next) = thresholdPeak;
            Dots.ItSum(next) = sum(thresholdMap(contourIndex));
            Dots.Vox(next).RawBright = iPunctaMed(contourIndex);
            Dots.MeanBrightPuncta(next) = mean(iPunctaMed(contourIndex));
            Dots.MeanBrightNeurite(next) = mean(iNeuriteMed(contourIndex));
            Dots.Contrast(next) = single(thresholdPeak)/single(mean(iPunctaMed(contourIndex)));
        end
                             
    else
        thresholdPeak = min(thresholdPeak);
        cutOff = 0.5 * thresholdPeak;
        contourIndex = find(thresholdLabel==i & thresholdMap>=cutOff);
        [yContour xContour zContour] = ind2sub([ys xs zs],contourIndex);
        unassignedContour = [yContour xContour zContour];
        nContours = numel(contourIndex);
        distance = zeros(nContours, nPeaks);
        for j=1:nPeaks
            distance(:,j) = sqrt(...
                (unassignedContour(:,1) - yPeak(j)).^2 +...
                (unassignedContour(:,2) - xPeak(j)).^2 +...
                (unassignedContour(:,3) - zPeak(j)).^2);
        end
        [dummy, peakAllegiance] = min(distance,[],2);
        for k=1:nPeaks
            peak(k).contour = unassignedContour(peakAllegiance==k,:); %#ok<AGROW>
            peak(k).contourIndex = contourIndex(peakAllegiance==k); %#ok<AGROW>
        end
        for l=1:nPeaks
            if  (yPeak(l)-1 < BlockBuffer/2 && yStart > 1) ||...
                    (xPeak(l)-1 < BlockBuffer/2 && xStart > 1) ||...
                    (zPeak(l)-1 < BlockBuffer/2 && zStart > 1) ||...
                    (ys-yPeak(l) < BlockBuffer/2 && yEnd < Rys) ||...
                    (xs-xPeak(l) < BlockBuffer/2 && xEnd < Rxs) ||...
                    (zs-zPeak(l) < BlockBuffer/2 && zEnd < Rzs)
            else
                if ~exist('Dots','var')
                    lastEntry = 0;
                else
                    lastEntry = size(Dots.Pos,1); 
                end
                next = lastEntry+1;
                Dots.Pos(next,:) = [yPeak(l)+yStart-1, xPeak(l)+xStart-1,...
                    zPeak(l)+zStart-1];
                Dots.Vox(next).Pos = [peak(l).contour(:,1)+yStart-1,...
                    peak(l).contour(:,2)+xStart-1,...
                    peak(l).contour(:,3)+zStart-1];
                Dots.Vox(next).Ind = sub2ind([Rys Rxs Rzs],...
                    Dots.Vox(next).Pos(:,1), Dots.Vox(next).Pos(:,2),...
                    Dots.Vox(next).Pos(:,3));
                Dots.Vol(next) = numel(peak(l).contourIndex);
                Dots.ITMax(next) = thresholdPeak;
                Dots.ItSum(next) = sum(thresholdMap(peak(l).contourIndex));
                Dots.Vox(next).RawBright = iPunctaMed(peak(l).contourIndex);
                Dots.MeanBrightPuncta(next) = mean(iPunctaMed(peak(l).contourIndex));
                Dots.MeanBrightNeurite(next) = mean(iNeuriteMed(peak(l).contourIndex));
                Dots.Contrast(next) = single(thresholdPeak)/single(mean(iPunctaMed(contourIndex)));
            end
        end
    end
end
        end
    end
end
Dots.ImSize = [Rys Rxs Rzs];
Dots.Num = size(Dots.Pos,1); 
save([TPN 'data/Dots.mat'],'Dots')
toc       
    

