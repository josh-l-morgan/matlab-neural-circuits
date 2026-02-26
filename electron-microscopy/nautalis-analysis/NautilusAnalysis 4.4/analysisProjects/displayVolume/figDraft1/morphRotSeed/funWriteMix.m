function[Imix] = funWriteMix(i,mix,iDegs,mixDeg,dirOrd,iNams,Imix,crop)

%     mixTarg = (find(degID==d,1));

Imix = zeros(crop(1,2)-crop(1,1)+1,crop(2,2)-crop(2,1)+1,3);

for m = 1:length(mix)
    if mix{m}(i)>0
        targ = find(iDegs{m}==mixDeg(i),1);
        if ~isempty(targ)
            
            for n = 0:10
                pass = 1
                try
                    I = double(imread([dirOrd{m} iNams{m}(targ+n).name]));
                catch err
                    pass = 0;
                end
                if pass, break,end
            end
            
            I = I(crop(1,1):crop(1,2),crop(2,1):crop(2,2),:);
            
            if exist([dirOrd{m} 'Label.png'],'file') & (mix{m}(i) == 1)
               Ilab = double(imread([dirOrd{m} 'Label.png']) );
                [ly lx lz] = size(Ilab);
                [iy ix iz] = size(Imix);
                Imix(1:min(ly,iy), 1:min(lx,ix),1: min(lz,iz)) = ...
                    Imix(1:min(ly,iy), 1:min(lx,ix),1: min(lz,iz)) + ...
                    Ilab(1:min(ly,iy), 1:min(lx,ix),1: min(lz,iz)) ;
            end
            
            Imix = Imix + I * mix{m}(i);
        end
        
    end
    
end
