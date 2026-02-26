

SPN = 'L:\joshm\LGNs1\PSC_alignments\S32intermediate\'
TPN = 'L:\joshm\LGNs1\PSC_alignments\S32intermediateSingleList\'

if ~exist(TPN,'dir')
    mkdir(TPN)
end


dSPN = dir(SPN); dSPN = dSPN(3:end);

for i = 1:length(dSPN)
    if dSPN(i).isdir
        nam = dSPN(i).name;
        disp(nam)
        dI = dir([SPN nam]); dI = dI(3:end);
        for f = 1: length(dI)
            tifName = dI(f).name;
            if sum(regexp(tifName,'.tif'))
                newName = [nam '_' tifName];
                if ~exist([TPN newName],'file')
                    for r = 1:3 % three tries at copy
                        pass = 1;
                        try  %try copying
                            copyfile([SPN nam '\' tifName],[TPN newName]);
                        catch err
                            err
                            pass = 0;
                        end
                        if pass  %copied ok
                            break
                        end
                    end
                    
                end
                
            end
        end
        
        
    end
    
end