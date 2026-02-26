%%Turn J OLE montage into fiji readable

TPN = GetMyDir;

dTPN = dir(TPN);
dTPN = dTPN(3:end);

iNam = [];
for i = 1:length(dTPN)
    nam = dTPN(i).name;
    if length(nam)>4
        if strcmp(nam(end-3:end), '.tif' )
            iNam{size(iNam,1)+1,1} = nam;
        end
    end
end

newTPN = [TPN(1:end-1) '_renamed\'];
if ~exist(newTPN)
    mkdir(newTPN)
end


%% find indicies
for i = 1:length(iNam)
   clear pos1
    sprintf('copying %d of %d',i,length(iNam))
   nam = iNam{i};
   vNam = nam;
   dNam = double(nam);
   dNums = (nam>=48) & (nam<=57);
   vNam(dNums) = 'V';
   pos1 = strfind(vNam,'_VVV_VVVxVVV');
   
   if ~isempty(pos1)
      tz = nam(pos1+1:pos1+3); 
      ty = nam(pos1+5:pos1+7);
      tx = nam(pos1+9:pos1+11);
      
      newNam = ['Tile_Z' tz '_Y' ty '_X' tx '.tif']
      pause
     % copyfile([TPN nam],[newTPN newNam]);
      
      
   end
    
    
end









