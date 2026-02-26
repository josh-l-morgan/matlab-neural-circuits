clear all
DPN = GetMyDir;
dDPN = dir(DPN); dDPN = dDPN(3:end);
sourceFolder = {}; maskFolder = {};
for d = 1:length(dDPN)
    if exist([DPN dDPN(d).name '\sourceFolder.mat']) 
        load([DPN dDPN(d).name '\sourceFolder.mat']) 
        sourceFolder{length(sourceFolder)+1} = getBip;
        maskFolder{length(maskFolder)+1} = dDPN(d).name;
    end
end
celldat = sourceFolder; %works relative to quickBip.m
%% Sort Data

[num txt dat] = xlsread('C:\Users\joshm\Documents\myData\FixedCon.xls','TomA');
%File	experiment date	Age	Retina	RGC	Bipolar	ID	description	Masked	to CB	CB rad	to center	center rad	bip class	rgc class	rgc depth	syn 	skeleton	use	notes											
colStrings = dat(1,:);
colFile = find(strcmp('File',colStrings));
colRet = find(strcmp('Retina',colStrings));
colRGC = find(strcmp('RGC',colStrings));
colBipolar = find(strcmp('Bipolar',colStrings));
colID = find(strcmp('ID',colStrings));

colStrings = dat(1,:);
colFile = find(strcmp('File',colStrings));
colRet = find(strcmp('Retina',colStrings));
colAge = find(strcmp('Age',colStrings));
colDistance = find(strcmp('distance',colStrings));
colBipClass = find(strcmp('bip class',colStrings));
colRgcClass = find(strcmp('rgc class',colStrings));
colSyn = find(strcmp('syn ',colStrings));
colNotes = find(strcmp('notes',colStrings));
colCenterRad = find(strcmp('center rad',colStrings));
colToCenter = find(strcmp('to center',colStrings));
colNumVox = find(strcmp('numVox',colStrings));
colNumTer = find(strcmp('numTer',colStrings));
colPolyArea = find(strcmp('polyArea',colStrings));
colVoxAppo = find(strcmp('voxAppo',colStrings));
colNumAppo = find(strcmp('numAppo',colStrings));

%% find the position of the cell in the data sheet using the file name
cellPos = zeros(length(celldat),1);
colList = [colFile colRet colRGC colBipolar colID];
for cD = 1: length(celldat) %find position of cell in data sheet
    clear namF
    %%break up name
    nam= celldat{cD}
    %nam = nam(1:end-1);
    slashes = find(nam=='\');
    L = length(slashes);    
    namF{4} = nam(slashes(L-2)+1:slashes(L-1)-1);
    namF{3} = nam(slashes(L-3)+1:slashes(L-2)-1);
    namF{2} = nam(slashes(L-4)+1:slashes(L-3)-1);
    namF{1} = lower(nam(slashes(L-5)+1:slashes(L-4)-1));
    namID = nam(slashes(L-1)+1:length(nam)-1);
    dash = find(namID =='_',1);
    if ~isempty(dash)
        namF{5} = lower(char(namID(dash+1:length(namID))));
    end
    for nf = 2:4
        n = lower(namF{nf});
        dash = find(n == '_',1);
        if isempty(dash)
            namF{nf} = n;
        else
            namF{nf} = n(1:find(n=='_',1)-1);
        end
    end
    
    
    %%search for name
    level = 1; i = 1;
    namF = lower(namF);
    memNam = [];
    while i < size(dat,1)     
        compTo = dat{i,colList(level)};
        compTo = lower(num2str(compTo));

        % i, namF{level}, compTo,  pause

        if strcmp(compTo,namF{level})
            level = level+1;
            memNam = [memNam '\' compTo];
            if level>length(namF)
                cellPos(cD) = i;
                break
            end
        else
            i = i + 1;
        end
    end
    memNam
    if ~cellPos(cD),pause,end
end

cP = 1:max(cellPos) * 0;
cP(cellPos) = 1:length(cellPos);
cP = cP';

tooMany = find((hist(cellPos,1:max(cellPos))>1))
tooFew = find((hist(cellPos,1:max(cellPos))<1))

find(cellPos == tooMany(2))

%% Read Dat
curFile = 'none';    curAge = 21;    curRet = 'none';    curRGC = 'none';
curRgcClass = 'none';
type6 = []; type7 = []; type8 = []; typeRb = []; typeU = []; typeShaft = [];
knownSyn = []; GU = [];
for i = 1:22 G{i} = []; end
for i = 1: size(dat,1)  %read data from sheet

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
    if ~isnan(dat{i,colRgcClass})
        curRgcClass = dat{i,colRgcClass};
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
            case 'shaft'
                typeShaft = [typeShaft ; i];
            otherwise
                typeU = [typeU ; i];
        end
    end

    if strcmp(curRgcClass(1), 'G');
        rc = str2num(curRgcClass(2:length(curRgcClass)));
        G{rc} = [G{rc} i];
    else
        GU = [GU i];
    end


    if ~strcmp(class(dat{i,colSyn}),'char')   & ~isnan(dat{i,colSyn})
        knownSyn = [knownSyn; i];
    end

    Files{i} = curFile;
    Ages(i) = curAge;
    Rets{i} = curRet;
    masked(i) = strcmp(class(dat{i,colPolyArea}),'double');
    appoed(i) = isnumeric(dat{i,colVoxAppo}) & sum(dat{i,colVoxAppo}>0) % strcmp(class(dat{i,colVoxAppo}),'double')
end

myG = [G{1} G{2} G{6} G{10}];
idTypes = {'type6', 'type7','typeRb', 'type8','typeShaft'};
col = {'r','g,','b','c','m'};
allAges = unique(Ages(~isnan(Ages)));
myMarkers = { 'o','x', '+', 's', 'd', 'p', 'h','v', '>','<'};

%% Translate xcell sheet data into typeID and ageID for cP
cAges = Ages(cellPos);
for i = 1:length(idTypes)
    typeID{i} = find(ismember(cellPos, eval(idTypes{i})));
end
ageID{1} = find(cAges == 9);
ageID{2} = find(cAges == 21);

%% read cP data
tips = zeros(length(maskFolder),1);
totLength = tips; polyArea = tips; depthArbor = tips; volume = tips;
branchPoints = tips; bulbosity = tips;
for i = 1:length(maskFolder)
    load([DPN maskFolder{i} '\data\cellProps.mat']);
    bulbosity(i,1) = cP.bulbosity;
    totLength(i,1) = cP.totLength;
    polyArea(i,1) = cP.polyArea;
    depthArbor(i,1) = length(cP.depthArbor);
    volume(i,1) = sum(cP.myProps.vol);
    links = cP.wire.links;
    hLinks = hist(links(:),min(links(:)):1:max(links(:)));
    tips(i,1) = sum(hLinks==1);
    branchPoints(i,1) = sum(hLinks>2);
end

%% Plot data
subplot(1,1,1)
mCol = { 'r' 'm'; 'b' 'c'; 'g' 'y'};
for i = 1:3
    for a = 2
        u = intersect(ageID{a}, typeID{i});
%         scatter3(polyArea(u),volume(u),bulbosity(u),'MarkerEdgeColor',mCol{i,a})
        
        subplot(3,2,1)
        scatter(tips(u),totLength(u),'MarkerEdgeColor',mCol{i,a},'Marker','.')
                hold on
        subplot(3,2,3)
        scatter(polyArea(u),depthArbor(u),'MarkerEdgeColor',mCol{i,a},'Marker','.')
        hold on
        subplot(3,2,5)
        scatter(bulbosity(u),volume(u),'MarkerEdgeColor',mCol{i,a},'Marker','.')

        hold on
    
    end
end
    hold off
pause(1)
cDat = [tips; totLength ; polyArea ; depthArbor;...
    bulbosity; volume];
useDat = [ typeID{1}; typeID{2}; typeID{3}];
useDat = intersect(useDat,ageID{2});
[IDX,C,sumd,D] = kmeans(cDat(useDat,:),3);

for i = 1:3
    for a = 2
       u = find(ismember(useDat, intersect(ageID{a}, typeID{i})));
       subplot(3,2,i*2)
       hist(IDX(u));
       hold on
    end
end
hold off


%% Check Sorting
colormap gray(256)

%for t = 1:length(typeID);
    t = 2
    tids = typeID{t};
    for ti = 1:length(tids)
    load([DPN maskFolder{tids(ti)} '\data\cellProps.mat']);
    wI = cP.wI;
    
    subplot(1,2,1);
    wSum = sum(wI>1,3);
    gam = (0:max(wSum(:))+1).^.5;
    gam = gam * 255/max(gam);
    image(gam(wSum+1)),
    
    subplot(1,2,2)
    wSum = squeeze(sum(wI>1,2))';
    gam = (0:max(wSum(:))+1).^.5;
    gam = gam * 255/max(gam);
    image(gam(wSum+1)),
    pause
    
end
%end


%%
emptyStrat = 0;
for t = 2
   useType = typeID{t};
   for i = 1: length(useType)
      strats = sortedStrat{useType(i)}; 
        if size(strats,2)>2
            
        plot(strats(:,1),'r'), hold on
        plot(strats(:,2),'b')
        plot(strats(:,3),'g')
        hold off 
            pause
        else
            emptyStrat = emptyStrat +1
        end
   
end

end