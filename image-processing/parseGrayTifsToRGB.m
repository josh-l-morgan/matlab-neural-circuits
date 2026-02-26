
%SPN = GetMyDir
SPN = 'F:\hxH\60x\60x_at_target_01.oif.files\'
TPN = [SPN(1:end-1) '_rgbTweak\'];
if ~exist(TPN,'dir'), mkdir(TPN);end


inams = dir([SPN '*.tif']);

clear C Z N
for i = 1:length(inams)
    nam = inams(i).name;
    cPos = regexp(nam,'C');
    zPos = regexp(nam,'Z');
    C(i) = str2num(nam(cPos(end)+1:cPos(end)+3));
    Z(i) = str2num(nam(zPos(end)+1:zPos(end)+3));
    N{i} = nam(1:max(1,cPos(end)-1));
end

%%
sectionList = unique(Z);
secNum = length(sectionList);
for i = 1:1:length(sectionList);
    i
    newName = sprintf('%s_%04.0f.tif',N{i},i);
    
            I = zeros(1,1,3,'double');

    for c = 1:3
        
        targ = find((C == c) & (Z==i),1);
        if exist([SPN inams(targ).name],'file')
            Ic =  imread([SPN inams(targ).name]);
            if 1%strcmp(class(Ic),'uint16');
                Ic = Ic/(16);
            end
            Ic = double(Ic);
            
            Ic = medfilt2(Ic,[3 3]);
            I(1:size(Ic,1),1:size(Ic,2),c) = Ic;
        end
    end
    
    
    %% tweak channel contrast
    It = I;
%     It(:,:,1) = (It(:,:,1)-10) * 2;
%     It(:,:,2) = (It(:,:,2)-10) * 2;
%     It(:,:,3) = (It(:,:,3)-20) * 10;
%     
    %% linear scale with depth
%     It = It * (1+((secNum-i)/500)*2.5);
    
    
    It = uint8(It);
    image(It),pause(.01)
    imwrite(It,[TPN newName]);

    
end



