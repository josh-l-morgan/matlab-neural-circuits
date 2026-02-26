function[typeLab typeKey] = getTypes(ids)

typeLab = [];
for i = 1:length(ids)
    nam = ids{i,2}
    ID = ids{i,1}
    if ID>0;
        if sum(regexp(lower(nam),'background'))
            typeLab(ID) = 0;
        elseif sum(regexp(lower(nam),'rgc'))
            typeLab(ID) = 1;
        elseif sum(regexp(lower(nam),'lin'))
            typeLab(ID) = 2;
        elseif sum(regexp(lower(nam),'tcp'))
            typeLab(ID) = 3;
        elseif sum(regexp(lower(nam),'gli'))
            typeLab(ID) = 4;
        elseif  sum(regexp(lower(nam),'segment'))
            typeLab(ID) = 5;
        elseif  sum(regexp(lower(nam),'mye'))
            typeLab(ID) = 6;
        elseif  sum(regexp(lower(nam),'xun'))
            typeLab(ID) = 7;
        else
            typeLab(ID) = 8;
        end
    end
    
end


typeKey = {'RGC';'Inhibitory'; 'thalamocortical'; 'glia'; ...
    'segment'; 'myelenated'; 'unknown';    'other'};
