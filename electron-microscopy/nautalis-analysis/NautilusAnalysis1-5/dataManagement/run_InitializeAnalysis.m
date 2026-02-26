
%%Select the directory to be analysed and add the nautalis function
%%directory and subdirectories to the matlab path



currentDirectory = pwd;
analysisFunctionsDirectory = [currentDirectory '\nautalisFunctions\'];
addpath(genpath(analysisFunctionsDirectory))

analysisFunctionsDirectory = [currentDirectory '\analysisProjects\'];
addpath(genpath(analysisFunctionsDirectory))

%MPN = GetMyDir;
