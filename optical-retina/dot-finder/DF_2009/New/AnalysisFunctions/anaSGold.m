function[] = anaSG(TPN)

if exist([TPN 'Settings.mat'])
        load([TPN 'Settings.mat'])
else
    Settings = [];
end

load([TPN 'Dots.mat'])

subplot(1,1,1)

if ~exist([TPN 'find']), mkdir([TPN 'find']); end
if ~exist([TPN 'images']), mkdir([TPN 'images']); end



%% Insure each Voxel corresponds to only one Dot
%%(Brutally removes all but brightest of overlapping dots)
%%Will cause problems if not deleted by rerun of dotfinder

IDs = 1:Dots.Num;
if exist([TPN 'data/Singles.mat'])
    load([TPN 'data/Singles.mat'])
        if length(Singles) ~= Dots.Num
            clear Singles
        else
            IDs = find(Singles);
        end
end

if ~exist('Singles')
    'Removing redundant dots...'
    ListInds = [];
    for i = 1: length(IDs)  % Find all indicies
        ListInds = [ListInds ; Dots.Vox(IDs(i)).Ind  ones(length(Dots.Vox(IDs(i)).Ind),1) * IDs(i)];
    end
    [uniqueInds uniqueFirst] = unique(ListInds(:,1),'first');
    [uniqueInds uniqueLast] = unique(ListInds(:,1),'last');
    repInds=uniqueInds((uniqueLast-uniqueFirst)>0);   %Find any repeated indicies

    %%remove dots that are not the brightest
    for i = 1: length(repInds)
        multiDot =[ListInds(find(ListInds(:,1)==repInds(i)),2)];
        howBright = Dots.ItSum(multiDot);
        brightest = multiDot(find(howBright == max(howBright),1));
        notBrightest = setdiff(multiDot,brightest);
        IDs = setdiff(IDs,notBrightest);
    end

    Singles = zeros(1,Dots.Num);
    Singles(IDs) = 1;  %Create logic matrix that excludes redundant dots
    save([TPN 'data/Singles.mat'],'Singles')
   
end


%% Define variable Defaults

%%Default dot criteria (first pass)
cP.distToMask = 1; %distance of dot peak to dend mask
cP.cut = 2; %ID from mask in which to remove dots
cP.compact = 10; % select for compact volumes
cP.contrast = 0.05; % measure of dot contrast (iterative threshold passes/ mean brightness)
cP.contrastVol = 1;  % Threshold for (ITMax*10)./(Dots.Vol.^(2/3));
cP.deltaFOf = 0.3; % Delta F over F (predicted for red) for puncta
cP.deltaF = 2; % Delta F (fluorecence - predicted fluorescence)
cP.vol = 3;  % Minimum volume to pass
cP.iTMax = 3; % Minimum times a dot must pass iterative threshold
cP.iTSum = 50; % Minimum sum of volume of each threshold pass
cP.deltaFOfTop = 0.3; % Minimum average delta F over f for the brightest %50 of voxels of a dot
cP.meanBright = 0; % Minimum mean brightness of puncta
cP.reset = 0; % Resets Criteria to defaults

defaultP.cP=cP; %Record as defaults

if isfield(Settings,'FirstThresholds')  % Check if FirstThresholds has been created
    cP = Settings.FirstThresholds; % Load FirstThresholds
else
    cP = defaultP.cP; % Use Default Thresholds
end

%% Get User input
% prompt = {'Volume in voxels','Maximum number of threshold passes',...
%     'Volume integrated threshold passes', 'Delta F/F (top half)','Reset Default'};
title = 'First Pass Thresholds: enter minimum values';
cP = getVars(cP,title);  % Let user define first criteria
if cP.reset, cP = defaultP.cP; end % If reset was pressed, reset values to defalut

%%Save Criteria
Settings.FirstThresholds = cP;
save([TPN 'Settings.mat'],'Settings')


%% Collect dot properties

clear dP
% dP.distToMask = Dots.DistToMask; faces the wrong direction
% dP.cut = Dots.Cut; %ID from mask in which to remove dots
% dP.compact = Dots.Compact; % select for compact volumes
dP.deltaFOf = Dots.DFOf; % Delta F over F (predicted for red) for puncta
dP.deltaF = Dots.DF; % Delta F (fluorecence - predicted fluorescence)
dP.vol = Dots.Vol;  % Minimum volume to pass
dP.iTMax = Dots.ITMax; % Minimum times a dot must pass iterative threshold
dP.iTSum = Dots.ItSum; % Minimum sum of volume of each threshold pass
dP.deltaFOfTop = Dots.DFOfTopHalf; % Minimum average delta F over f for the brightest %50 of voxels of a dot
dP.meanBright = Dots.MeanBright;

gBG=Dots.Im.GreenBackGround;  % Get the Background value of the Green Channel
mBG=dP.meanBright-gBG; % Get the Difference betteen dot brightness and background
dP.contrast=double(dP.iTMax)./dP.meanBright;% Define Contrast as number of times passing threshold divided by mean brightness of puncta  or %max(1,(mB-(It*2)-gBG));
dP.contrastVol=double(dP.iTMax*10)./(dP.vol.^(2/3)); % Scale Contrast according to volume of puncta, or   %Contrast=It./max(mBG,1);

dist2CB = Dots.Dist2CB;
if sum(isnan(dist2CB)), dist2CB = Dots.Pos(:,3); 'Using Z depth instead of CB', end

%% Applly Criteria to Puncta

runFields = fieldnames(dP); %threshold all dot properties
pass = ones(1,Dots.Num); %set up vector to record threshold passes
for i = 1: length(runFields)
    if isfield(cP,runFields{i}) % check if there is a threshold for that property
        pass = pass &   dP.(runFields{i})>=cP.(runFields{i}); %execute threhold
        %[runFields{i} ' ' num2str(sum(pass))]
    else
        ['no theshold for  ' runFields{i}] %notify if no threshold exists
    end
end
%%Apply distances threshold
pass = pass & (Dots.DistToMask <= cP.distToMask); % apply minimum distance
pass = pass & ~(Dots.Cut == cP.cut); % remove Cut ID
pass = pass & Singles; % Exclude redundant dots

P=find(pass')'; %% creat list of passing puncta
SG.pass1=pass';
save([TPN 'find\SG.mat'],'SG')

%% Draw image
%%Create relevant image matricies
maxID=zeros(Dots.ImSize(1),Dots.ImSize(2),'uint16');
AllmaxID=zeros(Dots.ImSize(1),Dots.ImSize(2),'uint16');
maxPassed=zeros(Dots.ImSize(1),Dots.ImSize(2),'uint8');
YXsize=Dots.ImSize(1)*Dots.ImSize(2); %?
DisAmAll=maxPassed; % Plot all Dot Voxels
DisAm=maxPassed;    % Plot voxels from Dots that pass criteria

%create pic of passed
for i =1:Dots.Num
    
    DisAmAll(mod(Dots.Vox(i).Ind-1,YXsize)+1)=DisAmAll(mod(Dots.Vox(i).Ind-1,YXsize)+1)+1;
    AllmaxID(mod(Dots.Vox(i).Ind-1,YXsize)+1)=i;
end

for i = 1: size(P,2)
    maxID(mod(Dots.Vox(P(i)).Ind-1,YXsize)+1)=P(i);
    maxPassed(mod(Dots.Vox(P(i)).Ind-1,YXsize)+1)=200; % Assign number ( could be dot value)
    DisAm(mod(Dots.Vox(P(i)).Ind-1,YXsize)+1)=DisAm(mod(Dots.Vox(P(i)).Ind-1,YXsize)+1)+1;  %Stack labels
end

%%Assign different values for pixels corresponding to one vs more then one dot
OverLapped=DisAm;
OverLapped(DisAm==1)=150;
OverLapped(DisAm>1)=255;


%load Raw
if exist([TPN 'images\maxRaw.mat'])
    load([TPN 'images\maxRaw.mat'])
elseif exist([TPN 'images\RawMax.tif']);
    maxRaw=imread([TPN 'images\RawMax.tif']);
    save([TPN 'images\maxRaw.mat'],'maxRaw');
else
    load([TPN 'Post.mat'])
    load([TPN 'Dend.mat'])
    maxRaw=max(Dend,[],3);
    maxRaw(:,:,2)=max(Post,[],3);
    save([TPN 'images\maxRaw.mat'],'maxRaw')
end

%%combine and save and image Comparison
TestID=uint16(maxRaw)*2^8;
TestCrit=maxRaw;
TestID(:,:,3)=maxID;
FindEmpty=sum(sum(TestCrit,1),2)==0;
if sum(FindEmpty), Empty=find(FindEmpty,1); else Empty=1; end
TestCrit(:,:,3)=OverLapped;
imwrite(TestCrit,[TPN 'find\TestCrit.tif'],'Compression','none')
%imwrite(TestID,[TPN 'images\TestID.tif'],'Compression','none')
image(TestCrit),pause(.01)

%%  Check if the user information has been provided, terminate if not
if ~exist([TPN 'find\yes.tif']) | ~exist([TPN 'find\no.tif'])
    'Please provide yes.tif and no.tif for further analysis'
end

%% Use user labeled images to Identify good and bad dots

Miss=[]; %Define vector to list artifacts
if exist([TPN 'find\no.tif'])
    NoIM=max(imread([TPN 'find\no.tif']),[],3);
    [NoLab NoNum]=bwlabel(NoIM>0);
    %%Identify Artifacts to eliminate

    for i = 1:NoNum
        foundID = maxID(find(NoLab==i));
        foundID = foundID(foundID>0);
        if isempty(foundID)
            foundID =AllmaxID(find(NoLab==i));
        end
        Miss=[Miss foundID'];
    end
    Miss = Miss(Miss>0);
    Miss = unique(Miss); %count each dot once.
end

Hit = []; %define vector to list identified puncta
if exist([TPN 'find\yes.tif'])
    YesIM=max(imread([TPN 'find\yes.tif']),[],3);
    [YesLab YesNum]=bwlabel(YesIM>0);
    for i = 1:YesNum
        FindYes=find(YesLab==i); %retreive IDs of dots to be eliminated
        if ~sum(DisAmAll(FindYes)>1) %if unambiguous
            Hit =[Hit AllmaxID(FindYes)'];
        end
    end
    Hit = Hit(Hit>0);
    Hit = unique(Hit);
end

%If no hits are defined, define everything that is not a miss as a hit
if isempty(Hit)
    Hit = unique(AllmaxID);
    Hit = setdiff(Hit,Miss);
    Hit = Hit(Hit>0);
end


SG.manual.Miss=Miss;
SG.manual.Hit=Hit;


%% Scale DFOfTopHalf by Dist2CB
polyD=polyfit(dist2CB(Hit),Dots.DFOfTopHalf(Hit)',1);
showRange=1:200;
showPoly=polyD(1)*showRange +polyD(2);
PredDelta=polyD(1)*dist2CB +polyD(2);

dP.deltaScale = (Dots.DFOfTopHalf * mean(Dots.DFOfTopHalf(Hit)))./PredDelta';

%% Combine Delta and Contrast
clear var
dP.deltaCon=(dP.deltaScale/var(dP.deltaScale(Hit))) .* (dP.contrast/var(dP.contrast(Hit)));

% %% Show measures
% showVars = {'deltaFOfTop' 'deltaScale' 'contrast' 'deltaCon'};
% for i = 1:length(showVars)
%     subplot(length(showVars)+1,1,i)
%     showY=dP.(showVars{i});
%     scatter(dist2CB,showY,'k','Marker','.','SizeData',20)
%     ylabel(showVars{i})
%     hold on
%     scatter(dist2CB(Hit),showY(Hit),'g')
%     scatter(dist2CB(Miss),showY(Miss),'r')
%     hold off
% 
% end
% subplot(length(showVars)+1,1,length(showVars)+1)
% pause(.01)


%% Pick best Thresholds

clear gP
gP.stepThresholds = 5;
gP.repSearch = 5;

% gP.compact = 0;
gP.deltaFOf = 0;
gP.deltaF = 0;
gP.vol = 0;
gP.iTMax = 1;
gP.iTSum = 0;
gP.deltaFOfTop = 1;
gP.meanBright = 0;
gP.contrast = 0;
gP.contrastVol = 0;
gP.deltaScale = 1;
gP.deltaCon = 0;

default.gP = gP; %Remember defaults

if isfield(Settings,'guideThresholds')  % Check if FirstThresholds has been created
    gP = Settings.guideThresholds; % Load FirstThresholds
end

pass2 = pass;  % Collect dots with tighter thresholds
if ~isempty(Miss)


    %%User defines which dot properties to use
    title = 'Which variables should change according to user labeling?';

    gP = getVars(gP,title);  % Let user define first criteria
    %if v.reset, v = deFault; end % If reset was pressed, reset values to defalut

    %%Save Criteria
    Settings.guideThresholds = gP;
    save([TPN 'Settings.mat'],'Settings')


    Step=gP.stepThresholds;
    clear Prop StartSearch StopSearch Interval Range

    
    
    %%Define variables to use
    checkProps=fieldnames(gP);
    useProps = {};
    for i =3:length(checkProps)
        if gP.(checkProps{i})
            useProps{length(useProps)+1} =checkProps{i};
        end
    end


    
    Pnum = length(useProps);
    if Pnum>0 % If using some properties to for user guidance
    ['Running ' num2str(gP.repSearch * Step^Pnum) ' combinations of thresholds.'] 
    clear TList
    for i = 1:length(useProps)
        Prop(i,:)=dP.(useProps{i});
    end

    %%%%%%%%%Replace steps with variance?

    for rep = 1:gP.repSearch %Search for best thresholds five times, narrowing search each time
        clear StartSearch StopSearch Range
        if rep==1
            StartSearch = double(min(Prop(:,Hit),[],2));
            StopSearch=double(max(Prop(:,Miss),[],2));
            Interval = double((StopSearch-StartSearch))/Step;
            clear Range
            for i = 1:size(Prop,1)
                Range(i,:)=StartSearch(i)-Interval(i):Interval(i):StopSearch(i)+Interval(i); %#ok<AGROW>
            end

        else
            StartSearch=min(meanThreshes'-Interval,[],2); %use previous interval as new range
            StopSearch=max(meanThreshes'+Interval,[],2);
            Interval=(StopSearch-StartSearch)/Step;
            for i = 1:size(Prop,1)
                Range(i,:)=StartSearch(i)-Interval(i):Interval(i):StopSearch(i)+Interval(i);
            end
        end

        %MemStartStop(:,:,rep) = [StartSearch StopSearch];

        %Create matix of threshold permuations
        siz=[];
        for i = 1:size(Range,1)
            siz=[siz size(Range,2)];
        end

        TList = zeros(prod(siz),Pnum,'uint8');
        ThreshList = double(TList)*0;
        for countIDs = 1:prod(siz) % run every position

            ndx = countIDs;
            if length(siz)<=Pnum,
                siz = [siz ones(1,Pnum-length(siz))];
            else
                siz = [siz(1:Pnum-1) prod(siz(Pnum:end))];
            end
            n = length(siz);
            k = [1 cumprod(siz(1:end-1))];
            for i = n:-1:1,
                vi = rem(ndx-1, k(i)) + 1;
                vj = (ndx - vi)/k(i) + 1;
                ThreshPosition(i,1) = vj;
                ndx = vi;
            end
            TList(countIDs,:) = squeeze(ThreshPosition);
        end


        clear FalseP FalseN ThreshList Eff
        for i = 1 : prod(siz)

            for runThresh = 1: Pnum
                currentThresh = Range(runThresh,TList(i,runThresh));
                FalseP(runThresh,:) = Prop(runThresh,Miss) >= currentThresh;
                FalseN(runThresh,:) = Prop(runThresh,Hit) < currentThresh;
                ThreshList(i,runThresh) = currentThresh;
            end

            FalsePsum=sum(sum(FalseP,1)==size(Prop,1)); % artifacts that pass all thresholds
            FalseNsum=sum(sum(FalseN,1)>0); % puncta that fail to pass one or more threshold
            Err=sum(FalsePsum + FalseNsum);
            Eff(i,:)=[FalsePsum FalseNsum Err];

        end  %Finish running threshold combinations


        BestThreshes=find(Eff(:,3)==min(Eff(:,3)));  % Find thresholds that made fewest mistakes
        Balance=Eff(BestThreshes,1)-Eff(BestThreshes,2); % Of best, what is the differnce betwen hits and misses
        BestBalance=find(abs(Balance)==min(abs(Balance))); % Of best, select those in which mistakes cancel out
        meanThreshes=mean(ThreshList(BestThreshes(BestBalance),:),1); % Pick Thresholds
        Cpos=BestThreshes(BestBalance(1));

        ErrorsMade=Eff(Cpos,3)  % Record erros for those thresholds

        %MemThresh.meanThreshes(rep,:)=meanThreshes
        %MemThresh.ErrorsMade(rep,:)=Eff(Cpos,:);

    end %repeat Threshold finding


    SG.Performance.AllKnown=size(Hit,2)+size(Miss,2);
    SG.Performance.KnownHits=size(Hit,2);
    SG.Performance.KnownMisses=size(Miss,2);
    SG.Performance.FalsePositives=Eff(Cpos,1);
    SG.Performance.FalseNegatives=Eff(Cpos,2);
    SG.Performance.Errors=Eff(Cpos,3);
    SG.Performance.ErrorRate=SG.Performance.Errors/SG.Performance.AllKnown;
    SG.Performance

%% Applly Criteria to Puncta
for i = 1 : length(useProps)
    pass2 = pass2 & (dP.(useProps{i}) >= meanThreshes(i));
end   

    end % if some properties were selected
end % if there are some misses





%run final pass adding manual Hits and Misses
passF=pass2;
if exist([TPN 'find\yes.tif'])
    passF(Hit)=1;
end
passF(Miss)=0;

SG.pass2=pass2';
SG.passF=passF';


%% Save second order properties
SG.SecOrder.DeltaScale=dP.deltaScale;
SG.SecOrder.Contrast=dP.contrast;
SG.SecOrder.ContrastV=dP.contrastVol;
SG.SecOrder.DeltaCon=dP.deltaCon;

save([TPN 'find\SG.mat'],'SG')


P=find(passF')'; %% list of passing puncta


%% Draw image
'Drawing images'

maxSum=zeros(Dots.ImSize(1),Dots.ImSize(2));
maxID1=maxSum; maxID2 = maxSum; maxID3 = maxSum;


for i = 1: size(P,2)

    maxID1(mod(Dots.Vox(P(i)).Ind-1,YXsize)+1)=rand*100+50;
    maxID2(mod(Dots.Vox(P(i)).Ind-1,YXsize)+1)=rand*100+50;
    maxID3(mod(Dots.Vox(P(i)).Ind-1,YXsize)+1)=rand*100+50;
    maxSum(mod(Dots.Vox(P(i)).Ind-1,YXsize)+1)=maxSum(mod(Dots.Vox(P(i)).Ind-1,YXsize)+1)+1;

end

maxPassed2=uint8((maxSum>0)*200);

OverLapped=maxSum*0;
OverLapped(DisAm>0)=60;
OverLapped(DisAm>1)=100;
OverLapped(maxSum>0)=OverLapped(maxSum>0)+150;

MaxC=maxID1+(maxSum>1)*1000;
MaxC(:,:,2)=maxID2+(maxSum>1)*1000;
MaxC(:,:,3)=maxID3+(maxSum>1)*1000;
MaxC=uint8(MaxC);
image(MaxC),pause(.01)
% image(max(fullID,[],3)),P(i),pause

%combine and save and image Comparison
SGCrit=maxRaw;
SGCrit(:,:,1)=SGCrit(:,:,1);
SGCrit(:,:,2)=SGCrit(:,:,2);
FindEmpty=sum(sum(TestCrit,1),2)==0;
if sum(FindEmpty), Empty=find(FindEmpty,1); else Empty=1; end
SGCrit(:,:,3)=OverLapped;
imwrite(SGCrit,[TPN 'find\SGCrit.tif'],'Compression','none')
imwrite(MaxC,[TPN 'find\MaxCids.tif'],'Compression','none')
%image(SGCrit), pause(.01)
%%

%% Plot Results
% hist(DeltaScale,0:.1:5)
% DeltaThresh
%
% hist(Contrast)
% ConThresh
%
% hist(DeltaCon)
% DConThresh
%
% hist(Dots.DFOfTopHalf)
%
%
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','b','EdgeColor','w')
% hold on
% hist(DeltaScale(DeltaScale<DeltaThresh),0:.1:5)
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','r','EdgeColor','w')
%
% hold off
%
%
% %
% pSGdeltaScale=DeltaScale>=DeltaThresh;
% pSGcontrast=Contrast>=ConThresh;
% pSGdcon=DeltaCon>=DConThresh;
% pSGflatDelta=Dots.DFOfTopHalf>=FlatDeltaThresh;

% end %find criteria if Yes and No are available




