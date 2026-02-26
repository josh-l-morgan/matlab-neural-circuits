
%%Select the directory to be analysed and add the nautalis function
%%directory and subdirectories to the matlab path

%C:\Users\joshm\Documents\MATLAB\jlm_Code\EM\Nautilus Analysis 1.0

disp('Initializing Nautilus Analysis')



currentDirectory = pwd;
analysisFunctionsDirectory = [currentDirectory '\nautilusFunctions\'];
addpath(genpath(analysisFunctionsDirectory))

analysisFunctionsDirectory = [currentDirectory '\analysisProjects\'];
addpath(genpath(analysisFunctionsDirectory))


analysisFunctionsDirectory = [currentDirectory '\describeTissue\'];
addpath(genpath(analysisFunctionsDirectory))

analysisFunctionsDirectory = [currentDirectory '\dataManagement\'];
addpath(genpath(analysisFunctionsDirectory))

analysisFunctionsDirectory = [currentDirectory '\otherPeoplesCode\'];
addpath(genpath(analysisFunctionsDirectory))

disp('Set directory containing MPN.mat, merge and analysis directories')
LPN = GetMyDir
addpath(LPN);
%MPN = GetMyDir;
