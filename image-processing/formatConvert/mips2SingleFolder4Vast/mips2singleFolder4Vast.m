


SPN = 'E:\affFull_fullRes_mipd\';
TPN = 'E:\affFull_fullRes_mipd_single\';
if ~exist(TPN,'dir'),mkdir(TPN);end

secDir = dir(SPN);
secFold = {secDir([secDir.isdir]>0).name};


for s = 1:length(secFold)
    sec = str2num(secFold{s})
    if ~isempty(sec)
        mipDir = [SPN secFold{s} '\0\'];
        
        if exist(mipDir,'dir')
            
            iDir = dir([mipDir '*.png']);
            iName = {iDir.name};
            
            
            for i = 1:length(iName)
                newName = sprintf('%s%d_%s',TPN,sec,iName{i});
                if ~exist(newName,'file')
                    copyfile([mipDir iName{i}],newName);
                end
            end
            
        end
    end
    
    
end