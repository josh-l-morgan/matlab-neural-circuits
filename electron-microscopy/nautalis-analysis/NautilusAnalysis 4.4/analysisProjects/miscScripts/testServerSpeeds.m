
MPN1 = 'X:\Active\joshm\LGNs1\mergeSeg_mat\' % ris jlmorgan
MPN2 = 'Y:\Active\Data\AT\Morgan Lab\joshm\LGNs1\mergeSeg_mat\' % ris fitzp
MPN3 = 'Z:\joshm\LGNs1\mergeSeg_mat\' % morganlab
clear risJLMorgan risFITZP morganlabTime
for i = 1:100
    
    tic
    load([MPN1 'dsObj.mat'])
    risJLMorgan(i) = toc
    
    tic
    load([MPN2 'dsObj.mat'])
    risFITZP(i) = toc
    
    tic
    load([MPN3 'dsObj.mat'])
    morganlabTime(i) = toc
    
    
    disp(sprintf('risJLM = %0.1f, risFITZP = %0.1f, morganlab = %0.1f',...
        risJLMorgan(end), risFITZP(end), morganlabTime(end)))
    
end