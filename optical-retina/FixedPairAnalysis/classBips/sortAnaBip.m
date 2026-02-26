clear all
%%Get directory with image folders
disp('getting folders')
subFolders = findFolders;
disp('got folders')
%% Get root
rootD = subFolders{1};
%mkdir([rootD '\imageIndex'])

%%Find image folders
ImageFolders={};
ImageDir = ImageFolders;
for i = 1:length(subFolders)
    Name=subFolders{i};
    slashes = find(Name == '\');
    lastFold = Name(slashes(length(slashes))+1:length(Name));
    if strcmp(lastFold(1:min(7,length(lastFold))),'bipMask')
            ImageFolders = [ImageFolders ; [Name '\']];
            ImageDir=[ImageDir ;{Name(1:length(Name)-length(lastFold))}];
    end
end

%% run all images
TotalImages=length(ImageFolders)
celldat = {};
for allI = 1:TotalImages

    TPN = ImageFolders{allI};
    cP = anaOne;
end



%% Sort Data

[num txt dat] = xlsread('\\Wongraid2\wonglab\JoshDaniel\FixedCon.xls','TomA');
%File	experiment date	Age	Retina	RGC	Bipolar	ID	description	Masked	to CB	CB rad	to center	center rad	bip class	rgc class	rgc depth	syn 	skeleton	use	notes											
colStrings = dat(1,:);
colFile = find(strcmp('File',colStrings));
colRet = find(strcmp('Retina',colStrings));
colRGC = find(strcmp('RGC',colStrings));
colBipolar = find(strcmp('Bipolar',colStrings));
colID = find(strcmp('ID',colStrings));
colList = [colFile colRet colRGC colBipolar colID];


cellPos = zeros(size(celldat,1),1);
for cD = 1: size(celldat,1)
    clear namF
    %% break up name
    nam= celldat{cD,1}
    slashes = find(nam=='\');
    L = length(slashes);    
    namF{4} = nam(slashes(L-2)+1:slashes(L-1)-1);
    namF{3} = nam(slashes(L-3)+1:slashes(L-2)-1);
    namF{2} = nam(slashes(L-4)+1:slashes(L-3)-1);
    namF{1} = lower(nam(slashes(L-5)+1:slashes(L-4)-1));
    namID = nam(slashes(L)+1:length(nam)-1);
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
    
    
    %% search for name
    level = 1; i = 1;
    while i < size(dat,1)     
        compTo = lower(dat{i,colList(level)});
        compTo = num2str(compTo);
            
%         namF{level}
%         compTo
%         pause
        if strcmp(compTo,namF{level})
            level = level+1;
            if level>length(namF)
                cellPos(cD) = i;
                break
            end
        else
            i = i + 1;
        end
    end
end

celldat{~cellPos,1}

%%Make list
for i = 1: length(cellPos)
    if cellPos(i)
        sortedDat(cellPos(i),:) = celldat(i,:);
    end
end







