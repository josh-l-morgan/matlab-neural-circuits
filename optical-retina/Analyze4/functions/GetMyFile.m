function[DFN DPN] = GetMyFile(showText)

if ~exist('showText')
    showText = 'Get File';
end

%% Get Folder 

if ~exist('.\temp'), mkdir('.\temp');end

if exist('.\temp\LastFile.mat')
     load(['.\temp\LastFile.mat']);
     if exist(LastFile.LastDPN)
        [DFN DPN]=uigetfile([LastFile.LastDPN LastFile.LastDFN],showText);
     else
        [DFN DPN] = uigetfile('*.*',showText);
     end
else
     [DFN DPN] = uigetfile('*.*',showText);
end

LastFile.LastDPN=DPN;
LastFile.LastDFN=DFN;
if LastFile.LastDPN>0
save('.\temp\LastFile.mat','LastFile')
end



