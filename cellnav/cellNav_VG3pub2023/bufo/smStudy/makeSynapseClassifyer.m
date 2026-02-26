function[s] = makeSynapseClassifyer(COI)

if ~exist('COI','var')
    SPN = [glob.datDir 'Analysis\Data\preproc\'];
    load([SPN 'COI.mat']);
end

%% define syn groups
cids = COI.cids;

clear s
s(1).cids = []; %Bipolar cell inputs (ribbons)
s(1).input = 1;
s(1).checkCids = 0;
s(1).synType = 2;
s(1).name = 'rib';

s(2).cids = []; %Conventional inputs
s(2).input = 1;
s(2).checkCids = 0;
s(2).synType = 1;
s(2).name = 'non-rib';

s(3).cids = cids(COI.isRGC); % Outputs to any RGC type
s(3).input = 0;
s(3).checkCids = 1;
s(3).synType = [];
s(3).name = 'RGCs';

s(4).cids = cids(COI.isAMC); % Outputs to any RGC type
s(4).input = 0;
s(4).checkCids = 1;
s(4).synType = [];
s(4).name = 'toAMCs';

%%Clean 1400 2005 2014 2021 3334?
s(5).cids = [1006 1022 1033 1052 1061 1065 1074 1075 1083 1084 1087 1088 ...
    1091 1096 1099 1100 1105 1111 1113 1117 1121 1126 1128 1137 1138 1142 ...
    1143 1153 1155 1160 1179 1195 1204 1208 1215 1228 1230 1232 1234 1238 ...
    1243 1258 1267 1271 1272 1279 1400 2005 2008 2009 2011 2013 2015 2018 ...
    2019 2020 2021 2030 2033 2100 2104 2106 2107 2110 3116 3306 3307 3308 ...
    3309 3310 3311 3323 3324 3325 3329 3333 3336 4013 4018 5314 6007 6012 ...
    6023 6024 6027 6035 6039 6070 6072 6077 6113 6127 6130 6137 6149]; % Outputs to any RGC type
s(5).input = 1;
s(5).checkCids = 1;
s(5).synType = [];
s(5).name = 'tracedBCs';

%%To unknown
s(6).cids = [0];
s(6).input = 0;
s(6).synType = [];
s(6).name = 'toUnk';

%%All outputs
s(6).cids = [];
s(6).input = 0;
s(6).checkCids = 0;
s(6).synType = [];
s(6).name = 'allOutputs';

%%All inputs
s(7).cids = [];
s(7).input = 1;
s(7).checkCids = 0;
s(7).synType = [];
s(7).name = 'allInputs';

%%All synapses
s(8).cids = [];
s(8).input = [];
s(8).checkCids = 0;
s(8).synType = [];
s(8).name = 'allSyn';

%%to Unk
s(9).cids = COI.unkCids;
s(9).input = 0;
s(9).checkCids = 1;
s(9).synType = [];
s(9).name = 'toUnk';

%%Bipolar cell sub types
L = length(s);
useBipNames = {'on' 'off' 'bc3a' 'bc3b' 'bc4' 'bc5i' 'bc5o' 'bc5t' 'xbc' 'bc6' 'bc7'};
c = 0;
for i = 1:length(COI.bpcGroupNames)
    name = COI.bpcGroupLabel{i};
    if sum(strcmp(useBipNames,name))
        c = c + 1;
        gCids = COI.bpcGroupCids{i};
        s(c+L).cids = gCids; % Outputs to any RGC type
        s(c+L).input = 1;
        s(c+L).checkCids = 1;
        s(c+L).synType = [];
        s(c+L).name = name;
    end
end

%%RGC cell sub types
L = length(s);
for i = 1: length(COI.rgcGroupCids)
    gCids = COI.rgcGroupCids{i};
    s(i+L).cids = gCids; % Outputs to any RGC type
    s(i+L).input = 0;
    s(i+L).checkCids = 1;
    s(i+L).synType = [];

%     name = [];
%     for b = 1:length(COI.rgcGroupNames{i})
%         name = [name  COI.rgcGroupNames{i}{b} ' '];
%     end
    s(i+L).name = COI.rgcGroupLabel{i};

end
