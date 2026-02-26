%% {'File','Age','Retina','RGC','Bipolar','distance','bip class','rgc class','syn ',[],'notes'%%
clear all

[num txt dat] = xlsread('\\Wongraid2\wonglab\JoshDaniel\FixedCon.xls','TomA');


%% find columns 

colStrings = dat(1,:);
colFile = find(strcmp('File',colStrings));
colRet = find(strcmp('Retina',colStrings));
colAge = find(strcmp('Age',colStrings));
colRGC = find(strcmp('Bipolar',colStrings));
colDistance = find(strcmp('distance',colStrings));
colBipClass = find(strcmp('bip class',colStrings));
colRgcClass = find(strcmp('rgc class',colStrings));
colSyn = find(strcmp('syn ',colStrings));
colNotes = find(strcmp('notes',colStrings));
colCenterRad = find(strcmp('center rad',colStrings));
colToCenter = find(strcmp('to center',colStrings));
colNumVox = find(strcmp('numVox',colStrings));
colNumArea = find(strcmp('numArea',colStrings));
colPolyArea = find(strcmp('polyArea',colStrings));

%% Read Dat+
curFile = 'none';    curAge = 21;    curRet = 'none';    curRGC = 'none';
type6 = []; type7 = []; type8 = []; typeRb = []; typeU = []; knownSyn = [];
for i = 1: size(dat,1)

    if ~isnan(dat{i,colFile})
        curFile = dat{i,colFile};
    end
    if ~isnan(dat{i,colAge}) & isnumeric(dat{i,colAge})
        curAge = dat{i,colAge};
    end
    if ~isnan(dat{i,colRet})
        curRet = dat{i,colRet};
    end
    if ~isnan(dat{i,colRGC})
        curRGC = dat{i,colRGC};
    end
    

    if isempty(dat{i,colBipClass})
        typeU = [typeU; i];
    else
    switch dat{i,colBipClass};
        case 6
            type6 = [type6 ; i];
        case 7
            type7 = [type7 ; i];
        case 8
            type8 = [type8 ; i];
        case 'rb'
            typeRb = [typeRb ; i];
        otherwise
            typeU = [typeU ; i];
    end
    end
    
    
    if ~strcmp(class(dat{i,colSyn}),'char')   & ~isnan(dat{i,colSyn})
        knownSyn = [knownSyn; i];
    end
    
    Files{i} = curFile;
    Ages(i) = curAge; 
    Rets{i} = curRet;  
    masked(i) = strcmp(class(dat{i,colPolyArea}),'double');
    
end


idTypes = {'type6', 'type7', 'type8','typeRb','typeU'};
col = {'r','g,','b','c','m'};
allAges = unique(Ages(~isnan(Ages)));


for a = 1:length(allAges)
for i = 1:length(idTypes)
assignin('base',idTypes{i}, intersect(eval(idTypes{i}),knownSyn)); %restric type list to known syn number examples
synNum = cell2mat(dat(intersect(eval(idTypes{i}),find(Ages==allAges(a))),colSyn));
% 
% subplot(2,1,1)
% scatter(ones(length(synNum),1)*i,synNum)
% xlim([-.5 length(idTypes)+ .5])

subplot(length(allAges),length(idTypes),i+ (a-1) * length(idTypes))

bins = (0:20);
hist(synNum,bins)
title(idTypes{i})
xlim([-0.5 20.5])
ylim([0 6])
synNums{i} = synNum;
meanNum(i) = mean(synNum);
varNum(i) = var(synNum);
scaledVar(i) = var(synNum)/mean(synNum);

end
end

%% Run Statistics
for i = 1: length(idTypes)
   idTypes{i}
   young = cell2mat(dat(intersect(eval(idTypes{i}),find(Ages<=12)),colSyn));
   old = cell2mat(dat(intersect(eval(idTypes{i}),find(Ages>12)),colSyn));
   ranksum(young,old)
        
end



%% Process ecentricities
%figure
toCent = cell2mat(dat(2:size(dat,1),colToCenter));
toCent = [0 ; toCent];
centRad = cell2mat(dat(2:size(dat,1),colCenterRad));
centRad = [0 ; centRad];
knownDist = find(~isnan(toCent));

useS= intersect(knownSyn,knownDist);
useNow = intersect(intersect(find(Ages == 9),useS),[type6; type7]);
%scatter(toCent(useNow),[dat{useNow,colSyn}],'r')  
plotBin(toCent(useNow),[dat{useNow,colSyn}],20,'r') 
hold on
 
useNow = intersect(intersect(find(Ages == 21),useS),[type6; type7]);;
plotBin(toCent(useNow),[dat{useNow,colSyn}],20,'b') 

% 
% %% fit to polynomial
% p = polyfit(toCent(useNow),[dat{useNow,colSyn}]',1);
% x = 0 : 1 : max(toCent);
% f = polyval(p,x);
% plot(x,f,'b')
% ylim([0,max([dat{useNow,colSyn}])]);
% hold off

%% Look at bipolar properties
for t = 1:4
    subplot(2,2,t)
    cType = idTypes{t};
    useNow = intersect(useS,find(Ages == 9));
    useNow = intersect(useNow,[eval(cType)]);
    useNow = intersect(useNow,find(masked));
    %useNow = intersect(useNow,useS([dat{useS,colSyn}]>0));
    X = [dat{useNow,colPolyArea}] * 0.0603^2 ;
    Y = [dat{useNow,colNumArea}];
    plotBin(X,Y,50,'r');
    title(cType)
    hold on

    useNow = intersect(useS,find(Ages == 21));
    useNow = intersect(useNow,[eval(cType)]);
    useNow = intersect(useNow,find(masked));
    %useNow = intersect(useNow,useS([dat{useS,colSyn}]>0));
    X = [dat{useNow,colPolyArea}] * 0.0603^2 ;
    Y = [dat{useNow,colNumArea}];
    plotBin(X,Y,50,'b');

    hold off
end



%% Compare synapses bip arbor
subplot(2,1,2)
useNow = intersect(useS,find(Ages == 9));
useNow = intersect(useNow,[type7]);
useNow = intersect(useNow,find(masked));
%useNow = intersect(useNow,useS([dat{useS,colSyn}]>0));
X = [dat{useNow,colPolyArea}] * 0.0603^2 ;
Y = [dat{useNow,colNumArea}];
plotBin(X,Y,50,'r');

hold on

useNow = intersect(useS,find(Ages == 21));
useNow = intersect(useNow,[type7]);
useNow = intersect(useNow,find(masked));
%useNow = intersect(useNow,useS([dat{useS,colSyn}]>0));
X = [dat{useNow,colPolyArea}] * 0.0603^2 ;
Y = [dat{useNow,colNumArea}];
plotBin(X,Y,50,'b');

hold off




%% Chart Cells
%%seperate by ages.  Chart bipolar type by 



