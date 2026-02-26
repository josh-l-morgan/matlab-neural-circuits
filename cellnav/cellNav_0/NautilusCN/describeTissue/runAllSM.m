load(['MPN.mat'])
load([WPN 'tis.mat'])

cids = tis.cells.cids;


for i = 120:length(cids)
    
    filename = sprintf('%sSMs\\sm_cid%d.mat',WPN,  cids(i));
    disp(sprintf('making %s',filename))
    
    if 1;%~exist(filename,'file')
        
        makeSM(cids(i))
    end
    
end