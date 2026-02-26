
load('MPN.mat')
load([MPN 'obI.mat'])
load([MPN 'dsObj.mat'])


allVox = [];
for i = 1:length(dsObj);
    allVox = cat(1,allVox,dsObj(i).subs);
end

%%

tracedVox = size(allVox,1);
potVox = 500 / .2 * 400 / .2 * 280 / .2;

tracedPC = tracedVox/potVox * 100