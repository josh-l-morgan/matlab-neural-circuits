function[Ifull] = showSubSums(vastSubs,sr,col);


%%

 allSubs = cat(1,vastSubs{:});
    fmin = double(min(allSubs,[],1));
    fsize = double(max(allSubs,[],1));
    clear allSubs
   
if exist('sr','var')
%     sr(:,1) = sr(:,1) - fmin(1);
%     sr(:,2) = sr(:,2) - fmin(2);
%     sr(:,3) = sr(:,3) - fmin(3);
sr = round(sr);

 fmin = sr(1,:);
 fsize = sr(2,:);
  fsize = fsize-fmin+1;
else
    sr = cat(1,[0 0 0],fsize*2);
    
    
 fmin = max(fmin,sr(1,:));
 fsize = min(fsize,sr(2,:));
 fsize = fsize-fmin+1;
    
end



 


obNum = length(vastSubs);
if ~exist('col','var')
if obNum >1
    colMap = hsv(256);
    col = colMap(ceil((1:obNum)*256/obNum),:);
    col = col(randperm(size(col,1)),:);
else
    col = [1 1 1];
end
end

%%

if ~exist('dim','var')
    dim = 1;
end



%%
if ~exist('viewType','var')
    viewType = 1;
end

if viewType == 1
    
    minInt = 20;
    contrastFactor = 10;
    %%
    
    for d = 1: 3;
        dim = d;
        
        if dim == 1
            dims = [3 2];
        elseif dim == 2
            dims = [3 1];
        elseif dim == 3
            dims = [1 2];
        end
        
        %fsize = [1700 1700 1300];
        I = zeros(fsize(dims));
        Ic = cat(3,I,I,I);
        
        uCellId =  1:obNum;
        
        for i = 1:length(uCellId)
            
            sub = double(vastSubs{i});
            pass = (sub(:,1)> sr(1,1) ) & (sub(:,1)<sr(2,1) ) &...
                (sub(:,2)> sr(1,2) ) & (sub(:,2)< sr(2,2) ) &...
                (sub(:,3)> sr(1,3) ) & (sub(:,3)< sr(2,3) );
            sub = sub(pass,:);
            
            
            if ~isempty(sub)
                inds = sub2ind(fsize(dims),sub(:,dims(1))-fmin(dims(1))+1,sub(:,dims(2))-fmin(dims(2))+1);
                uinds = unique(inds);
                
                if length(uinds)>1
                    hinds = hist(inds,uinds);
                else
                    hinds = length(inds);
                end
                
                I(uinds) = I(uinds) + hinds';
                % newHeight(uinds) = sub(:,dim);
                
            
            Ibackup = I;
            %%
            I = Ibackup;
            GammaValue = .5;
            
            I = (255 - minInt) * (I/max(I(:))).^ GammaValue;
%             pause(.01)
%             image(I)
            
            %%
            
            %I = I * contrastFactor;
            I((I>0)) = I(I>0) + minInt;
            I(I>255) = 255;
            I((I>0) & (I<minInt)) = minInt;
            
            
            Ic(:,:,1) = Ic(:,:,1) + I*col(i,1);
            Ic(:,:,2) = Ic(:,:,2) + I*col(i,2);
            Ic(:,:,3) = Ic(:,:,3) + I*col(i,3);
            
            end
%             image(uint8(Ic))
%             'color'
%             pause(.01)
            I = I * 0;
            %%
            
        end
        Ic = uint8(Ic);
        Iall{d} = Ic;
        Isize{d} = fsize(dims);
        
    end
    
    
end


%% full image
% 
% Ifull = zeros(fsize(1) + fsize(3), fsize(2) + fsize(3),3,'uint8');
% Ifull(1:fsize(1),1:fsize(2),:) = Iall{3};
% Ifull(fsize(1)+1:fsize(1)+fsize(3),1:fsize(2),:) = Iall{1};
% Ifull(1:fsize(1),fsize(2)+1:fsize(2)+fsize(3),:) = permute(Iall{2},[2 1 3]);
% 
% % image(Ifull)



%%
Ifull = zeros(fsize(1) + fsize(3), fsize(2) + fsize(1),3,'uint8');
Ifull(fsize(3) + 1:fsize(3) + fsize(1),1:fsize(2),:) = Iall{3};
Ifull(1:fsize(3),1:fsize(2),:) = Iall{1};
Ifull(1:fsize(3),fsize(2)+1:fsize(2)+fsize(1),:) = Iall{2};
Ifull = permute(Ifull,[2 1 3]);

Ifull = uint8(Ifull);
% image(Ifull)
%imwrite(Ifull,[MPN 'threeViews.png'])

