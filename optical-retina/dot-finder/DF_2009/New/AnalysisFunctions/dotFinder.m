function[]=dotFinder(TPN) %#ok<INUSL>
% 
% TPN = GetMyDir;
%large files will be broken into smaller blocks and then recombined


%% Get file names
tic

load([TPN 'Settings.mat'])
load([TPN 'Post.mat']);
[Rys Rxs Rzs] = size(Post);
image(max(Post,[],3)),pause(.01)

colormap gray(255) %standard grey colormap
if isdir([TPN 'temp'])==0, mkdir([TPN 'temp']); end %create directory to store steps
if isdir([TPN 'data'])==0, mkdir([TPN 'data']); end %create directory to store steps
if isdir([TPN 'pics'])==0, mkdir([TPN 'pics']); end %create directory to store steps
 


%% Enter Variables

v = Settings.dotfinder;
% 
% %%Image Variables
% % xyum=.103;
% % zum=.3;
% % % aspect=zum/xyum;% ratio of z to xy dimentions
% % channels=3;
% v.blockSize=200; %aproximate size of individual processing blocks
% v.blockBuffer=30; %amount of block to be cut off as edge buffer
% 
% %%Dot criteria
% step=2; %Sensitivity=grey value step of iterative threshold (2 was standard)
% MaxDot=7^3;  %Maximum Dot Volume (in pixels)= maximum dot size for iterative threshold (6^3 was standard)
% MinDot=5;  %Minimum Dot Volume (in pixels)= minimum dot size for iterative threshold   (10 was standard)
% %MinDot=3;  %Minimum Dot Volume (in pixels)= minimum dot size for iterative threshold   (10 was standard)
% % PunctaThreshold=1; %Minimum Number of Steps Passed = centroids less than puncta threshold are zeroed. 
% % EdgeOfPeak=.5; %Determine Dot Edge = ratio of edge brightness to peak brightness (.5 was standard)
% % minFilledVolume=3; %minimum number of pixels in final object contour (20 was standard)
% % RoundThreshold=50; %minimum roundness threshold (60 was standard)


%% Find Image Stats
IgStat = sort(Post(:));
%IgStat = IgStat(IgStat>0);
lowpoint=min(IgStat(IgStat>0)); %find top and bottom
highpoint=IgStat(length(IgStat)); 
%%Determine intensity value at which 'v.percentBackground' is excluded from search
Gmode = IgStat(round(v.percentBackground*length(IgStat)));
if Gmode<lowpoint,Gmode = lowpoint; end
% Gmode = mode(single(IgStat(IgStat>0)));


%% SubSample
if Rxs>v.blockSize
    NumBx=round(Rxs/v.blockSize);
    Bxc=fix(Rxs/NumBx); 
else
    Bxc=Rxs;
    NumBx=1;
end

if Rys>v.blockSize
    NumBy=round(Rys/v.blockSize);
    Byc=fix(Rys/NumBy); 
else
    Byc=Rys; 
    NumBy=1; 
end

if Rzs>v.blockSize
    NumBz=round(Rzs/v.blockSize);
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
yStart=Tystart-v.blockBuffer;
yStart(yStart<1)=1;
yEnd=Tyend+v.blockBuffer;
yEnd(yEnd>Rys)=Rys;
xStart=Txstart-v.blockBuffer;
xStart(xStart<1)=1;
xEnd=Txend+v.blockBuffer;
xEnd(xEnd>Rxs)=Rxs;
zStart=Tzstart-v.blockBuffer;
zStart(zStart<1)=1;
zEnd=Tzend+v.blockBuffer;
zEnd(zEnd>Rzs)=Rzs;

load([TPN 'Post.mat']);
Igm = single(Post(yStart:yEnd,xStart:xEnd,zStart:zEnd)); 
clear Post
[ys,xs,zs]=size(Igm);


%% FIND DOTS Green Channel%%
peakMap = zeros(ys,xs,zs,'uint8');  %set up matrix to map peaks
thresholdMap = zeros(ys,xs,zs,'uint8');   %set up matrix to sum passed thresholds
% indVect = 1:ys*xs*zs;
maxIntensity = uint8(max(Igm(:)));

for i = maxIntensity:-v.thresholdStep:Gmode
    %run thresholds through all relevant intensities
    clear Igl labels
    [Igl,labels] = bwlabeln(Igm>i,6);%label each area to check

    %reduce bitdepth if possible
    if labels<65536
        Igl=uint16(Igl);
    end
    if labels <= 1
        labels =2;
    end
    nPixel = hist(Igl(Igl>0), 1:labels);
    %run all lables
    for p=1:labels
        % Morphology Filter, Puncta size criteria
        pixelIndex = find(Igl==p);
        if nPixel(p) < v.maxDotSize && nPixel(p) > v.minDotSize
            %identify peak in labeled object
            if sum(peakMap(pixelIndex))== 0
                peakValue = max(Igm(pixelIndex));
                peakIndex = find(Igl==p & Igm==peakValue);
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
        if  (yPeak-1 <= v.blockBuffer/2 && yStart > 1) ||...
                (xPeak-1 <= v.blockBuffer/2 && xStart > 1) ||...
                (zPeak-1 <= v.blockBuffer/2 && zStart > 1) ||...
                (ys-yPeak < v.blockBuffer/2 && yEnd < Rys) ||...
                (xs-xPeak < v.blockBuffer/2 && xEnd < Rxs) ||...
                (zs-zPeak < v.blockBuffer/2 && zEnd < Rzs)
        else
            cutOff = v.peakCutoff * thresholdPeak;
            contourIndex = find(thresholdLabel==i & thresholdMap>=cutOff);
            [yContour xContour zContour] = ind2sub([ys xs zs],contourIndex);
            if ~exist('Dots','var')
                lastEntry = 0;
            else
                [lastEntry dummy] = size(Dots.Pos);  %#ok<NASGU>
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
            Dots.Vox(next).RawBright = Igm(contourIndex);
            Dots.MeanBright(next) = mean(Igm(contourIndex));
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
            if  (yPeak(l)-1 <= v.blockBuffer/2 && yStart > 1) ||...
                    (xPeak(l)-1 <= v.blockBuffer/2 && xStart > 1) ||...
                    (zPeak(l)-1 <= v.blockBuffer/2 && zStart > 1) ||...
                    (ys-yPeak(l) < v.blockBuffer/2 && yEnd < Rys) ||...
                    (xs-xPeak(l) < v.blockBuffer/2 && xEnd < Rxs) ||...
                    (zs-zPeak(l) < v.blockBuffer/2 && zEnd < Rzs)
            else
                if ~exist('Dots','var')
                    lastEntry = 0;
                else
                    [lastEntry dummy] = size(Dots.Pos); %#ok<NASGU>
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
                Dots.Vox(next).RawBright = Igm(peak(l).contourIndex);
                Dots.MeanBright(next) = mean(Igm(peak(l).contourIndex));
            end
        end
    end  % if peaks
end  %All labels
        end
    end
end
Dots.ImSize = [Rys Rxs Rzs];
Dots.Num = size(Dots.Pos,1); %#ok<NASGU> 
save([TPN 'Dots.mat'],'Dots')
toc       
    

