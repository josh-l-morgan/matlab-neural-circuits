%% Run auto seg on stack

TPN = GetMyDir;

inams = getPics(TPN);  %find all tifs

if ~exist([TPN 'labeled'])
    mkdir([TPN 'labeled']);
end

for i = 1:length(inams)
   sprintf('running plane %d of %d',i,length(inams))
   %if ~exist([TPN 'labeled\lab' inams{i}])
   I = imread([TPN inams{i}]);    
   %I = I(600:900,600:900);
    I = segEWA(I);
    imwrite(uint16(I),[TPN 'labeled\lab' inams{i} '.tif'],'Compression','none');
   %end
end