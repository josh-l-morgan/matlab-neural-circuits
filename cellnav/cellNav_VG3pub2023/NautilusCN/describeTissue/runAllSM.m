function[] = runAllSM()

global globSM



load(['MPN.mat'])
load([WPN 'tis.mat'])

cids = tis.cells.cids;


for i = 1:length(cids)
    str = sprintf('running cell %d (%d of %d)',cids(i),i,length(cids));
    disp(str)
    try
        set(globSM.handles.textOut,'String',str);
    end
    filename = sprintf('%sSMs\\sm_cid%d.mat',WPN,  cids(i));
    disp(sprintf('making %s',filename))
    
    if 1;%~exist(filename,'file')
        
        makeSM(cids(i))
    end
    
end