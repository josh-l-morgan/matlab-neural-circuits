%%Give each color channel of a OIF its own directory of images
%%To make it easier to import into VAST

%% Set variables

convertTo8bit = 1;
medianFilter = 2; %kernal size, 0 for none


%% Get file names
%SPN = uigetdir;
%SPN = 'I:\LxA_GadGFP_Albino\Mouse4\SectionB4\60x2xStack.oif.files';
SPN = 'D:\Ziyi Hu\Mouse4\SectionB4\60x2xStack.oif.files';
TPN = [SPN '_ChannelFolders'];
if ~exist(TPN,'dir'),mkdir(TPN);end

dSPN = dir([SPN '\*.tif']);
nams = {dSPN.name};


%% Parse file names
z = zeros(length(nams),1);
c = z;
for i = 1:length(nams)
    nam = nams{i};
    cPos = regexp(nam,'C');
    zPos = regexp(nam,'Z');
    suf = regexp(nam,'.tif');
    zstr = nam(zPos(1)+1:suf(1)-1);
    z(i) = str2num(zstr);
    cstr = nam(cPos(1)+1:zPos(1)-1);
    c(i) = str2num(cstr);
end

%% Make channel directories
chans = unique(c);
for i = 1:length(chans);
    chanDir{i} = [TPN '\' num2str(chans(i))];
    if ~exist(chanDir{i},'dir'),mkdir(chanDir{i});end
end


%% Copy files to channel directories

if convertTo8bit


    for i = 1:length(nams)
        disp(sprintf('moving and converting file %d of %d',i, length(nams)))
        oldFile = [SPN '\' nams{i}];
        newFile = [chanDir{c(i)} '\' nams{i}];
        I = imread(oldFile);
        I = double(I);
        if medianFilter
            I = medfilt2(I,[ medianFilter medianFilter]);
        end
        I = I * (2^8 / 2^12);
        I = uint8(I);
        imwrite(I,newFile);
    end

else

    for i = 1:length(nams)
        disp(sprintf('copying file %d of %d',i, length(nams)))
        oldFile = [SPN '\' nams{i}];
        newFile = [chanDir{c(i)} '\' nams{i}];
        copyfile(oldFile,newFile);
    end

end





