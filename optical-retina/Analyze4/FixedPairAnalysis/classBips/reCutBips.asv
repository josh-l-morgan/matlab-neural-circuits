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
bipDump = 'C:\Users\joshm\Documents\myData\bipolarData\bipDump\'
for allI = 1:TotalImages
    getBip = ImageFolders{allI};
    nam = num2str(allI);

% 
%     copyfile(getBip,[bipDump nam],'f');
%     save([bipDump nam '\sourceFolder.mat'],'getBip')

end
