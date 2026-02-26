%%
clear all

%CPN = 'F:\exY_R1Cerebellum\'
CPN = 'W:\'

dCPN = dir([CPN ]);

c = 0;
for i = 1:length(dCPN)
   
    WPN = [dCPN(i).name '\'];
    dWPN = dir([CPN WPN '*_Montage']);
    
    for w = 1:length(dWPN)
        
        SPN = [dWPN(w).name '\'];
        dSPN = dir([CPN WPN SPN 'Tile_*.tif']);
        
        for t = 1:length(dSPN)
           
            nam = dSPN(t).name;
            
            rp = regexp(nam,'_r');
            cp = regexp(nam,'-c');
            wp = regexp(nam,'_w');
            sp = regexp(nam,'_sec');
            dp = regexp(nam,'.tif');
            
            if isempty(regexp(nam(sp+4:dp-1),'_'))
                c = c+1
                nam
                'meow'
            path{c,1} = [CPN WPN SPN nam];
            row{c,1} = nam(rp+2:cp-1);
            col{c,1} = nam(cp+2:wp-1);
            waf(c) = str2num(nam(wp+2:sp-1));
            sec(c) = str2num(nam(sp+4:dp-1));
            else
                'bark'
            end 
        end 
        
    end
    
end

id = waf * 10000 + sec;
uid = unique(id);
[sortID idx] = sort(uid,'ascend');

for i = 1:length(id)
   z(i) = find(sortID== id(i)); 
end


%%
% 
% pathString = [];
% posString = [];
% fileID = fopen([CPN 'imagePath2.txt'],'w');
% 
% for i = 1:length(path)
%     fprintf(fileID,'%s\r\n',path{i});
% end
% % fclose(fileID)
% 
% 
% fileID = fopen([CPN 'imagePos2.txt'],'w');
% for i = 1:length(path)
%     fprintf(fileID,'%s %s %d \r\n', col{i}, row{i},z(i));
% end
% fclose(fileID)


fileID = fopen([CPN 'imagePathPos.txt'],'w');
for i = 1:length(path)
    fprintf(fileID,'%s %s %s %d \r\n', path{i},col{i}, row{i},z(i));
end
fclose(fileID)

%%
