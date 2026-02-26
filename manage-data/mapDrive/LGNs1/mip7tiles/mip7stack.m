SPN = '\\storage1.ris.wustl.edu\jlmorgan\Active\morganLab\DATA\LGN_Adult\LGNs1\Josh_LGN_aligned_v3\em\mip7\'
TPN = '\\storage1.ris.wustl.edu\jlmorgan\Active\morganLab\DATA\LGN_Adult\LGNs1\Josh_LGN_aligned_v3\mip7stack\'

dSPN = dir([SPN 'slice*'])
slices = {dSPN.name};
numSlices = length(slices);

for i = 1:numSlices
    disp(sprintf('Copying %d of %d',i,numSlices));
    dI = dir([SPN slices{i} '\*.png']);
    fileName = dI(1).name;
    if ~exist([TPN fileName],'file')
    copyfile([SPN slices{i} '\' fileName],[TPN fileName]);
    end
end
