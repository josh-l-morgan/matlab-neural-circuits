function[] = tt_saveImageFunction


global globTT



TPN = globTT.dir.saveFolder;
TFN = globTT.dir.saveName;
overWrite = globTT.save.overWrite;

id = globTT.save.tab;
if id
    %I = globTT.I.tab{id};
    useID = id;
    useMix = 0;
    alph = ones(6,1);
else
    %I = globTT.save.I;
    useID = globTT.imMix.useID;
    alph = globTT.imMix.alpha;
    useMix = 1;
end

doStack = globTT.save.doStack; %& (size(I,4)>1);



%% Save
if ~doStack %write single image
    slice = globTT.active.slice;
     NFN = TFN;
    if ~overWrite
        for f = 1:999
            if exist([TPN NFN '.png'],'file')
                NFN = sprintf('%s_v%03.0f',TFN,f);
            else
                break
            end
        end
    end
    %if id, I = applyIdContrast(I); end
    I = applyMixAlpha(slice,useID,alph);
    
    imwrite(uint8(I),[TPN NFN '.png'],'WriteMode','overwrite');
    
    
else  %write stack into folder
    NFN = TFN;
    if ~overWrite
        for f = 1:999
            if exist([TPN NFN],'dir')
                NFN = sprintf('%s_v%03.0f',TFN,f);
            else
                break
            end
        end
    end
    
    secNum = size(globTT.I.tab{1},4);
    
    if ~exist([TPN NFN],'dir'), mkdir([TPN NFN]); end
    for i = 1:secNum 
        
        planeName = sprintf('%s%s\\%04.f.png',TPN,NFN,i);
        disp(['Writing ' planeName])
        %Ip = I(:,:,:,i);
        %if id, Ip = applyIdContrast(Ip); end
        Ip = applyMixAlpha(i,useID,alph); 
        imwrite(uint8(Ip),planeName,'WriteMode','overwrite');
    end
    
end


end

function[I] = applyIdContrast(I,id)
global globTT

for c = 1:3
    F = I(:,:,c);
    F = F + globTT.twk.bright(id,c);
    F = F * globTT.twk.con(id,c);
    F(F<0) = 0;
    maxA = max(F(:));
    conA =  maxA/(maxA ^ globTT.twk.gamma(id,c));
    
    F = F.^globTT.twk.gamma(id,c) * conA; %gamma presserving max value.
    I(:,:,c) = F;
end
end

function[Im] = applyMixAlpha(slice,useID,alph);

global globTT

Im = zeros(size(globTT.I.tab{useID(1)},1),size(globTT.I.tab{useID(1)},2),3);
for i = 1:length(useID)
    id = useID(i);
    %I = eval(sprintf('globTT.I.%s(:,:,:,%d);',globTT.active.I,globTT.active.slice));
    if size(globTT.I.tab{id},4) >= slice
        I =  globTT.I.tab{id}(:,:,:,slice);
    else
        I =  globTT.I.tab{id}(:,:,:,1);
    end
    
    I = applyIdContrast(I,id);
    
    Im  = Im + I * alph(id);
end

end






