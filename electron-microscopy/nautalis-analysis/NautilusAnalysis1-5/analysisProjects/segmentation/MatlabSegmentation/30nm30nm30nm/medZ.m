
%matlabpool open 7


colormap(gray(256))
SPN = GetMyDir;
TPN = [SPN(1:end-1) '_medZ' '\'];
mkdir(TPN)

%% read names
dSPN = dir(SPN);
iNam = {};
for i  =   1:length(dSPN)
   nam = dSPN(i).name;
   if sum(regexp(nam,'.tif'))
      iNam{length(iNam)+1} = nam; 
   end    
end


%%

iBin = 1:16:length(iNam);

iSize = 800;

parfor i = 1:length(iBin)-1;
    sI = zeros(iSize,iSize,4,'single');

    p = 0;
    disp(sprintf('creating image %d of %d',i,length(iBin)))
    for s = iBin(i):iBin(i+1)-1;
        
        for r = 1:3
            try
                pass = 1;
                I = imread([SPN iNam{s}]);
            catch err
                pass = 0
            end
            if pass
                p = p+1;
                sI(:,:,p) = single(I);
                break
            end
        end
    end
        
        mI = uint8(median(sI(:,:,1:p),3));
        %image(mI),pause(.01)
        startName = iNam{iBin(i)};
        stopName = iNam{iBin(i+1)-1};
        newName = [startName(1:end-4) '-' stopName(5:end)];
        
        imwrite(mI,[TPN newName]);
     
end



