function [masks,seeds] = calciumRegion(image,seeds,Ths,SEs,W)
%extract the mask and final seeds for cells
%   Inputs:
%     image is a 2D array. 
%
%     seed is a 2D array. Each row shows a seed
%     coordination.
%
%     Ths is a vector contaning the initial thresholds for the region growing segmentation starting from the seed points. 
% 
%     SEs is a cell array containing the Structural Elements for the close morphological operations. 
%     W is the window size for the cells around the seeds.
%      
%  Outputs:
%        masks is a 3D cell array. the first dimention is for the Ths, the
%        second one is for the seeds, and the last one is for the SEs.

im=double(image);
ims=size(image);
if isempty(seeds)
    % an auto seeds detection algorithm comes here.
    seeds=round(size(image)/2);
end
if isempty(Ths)
    Ths=0.001;
end
if isempty(SEs)
    SEs{1}=[0 1 0; 1 1 1; 0 1 0];
end
if isempty(W)
    W=0;
end
th=0.01; % threshold for the second region growing algorithm. the auto calculation comes here.
si=size(seeds);
masks = cell(size(Ths,2),si(1,1),size(SEs,2)); 
for it=1:size(Ths,2)% for each threshold
    for i=1:si(1,1)% for each seed point
        %% calculate a window around the seed point
        if W>0
        x1=seeds(i,1)-W;x2=seeds(i,1)+W;
        y1=seeds(i,2)-W;y2=seeds(i,2)+W;
        if x1<1
            x1=1;
        end
        if x2>ims(1)
            x2=ims(1);
        end
        if y1<1
            y1=1;
        end
        if y2>ims(2)
            y2=ims(2);
        end
        im=image(x1:x2,y1:y2);
        end
        %%
        %% calculate the seed location in the window
        ss=seeds(i,:)-W;
        if ss(1)<1
            ss(1)=1;
        end
        if ss(2)<1
            ss(2)=1;
        end
        %%
        %% apply region growing algorithm in the window starting from the seed point ss and compute the first mask
        Jrg=regiongrowing(im,ss(1),ss(2),Ths(it));
        %% correct the first mask by applying the close morphological operation with each structural element
        for is=1:size(SEs)
                    Jse = imclose(Jrg,SEs{is});
                    Je = edge(Jse,'sobel');
                    [x,y]=ind2sub(size(im),find(Je));
                    
                    for k=1:size(x)
                        if Jse(x(k),y(k))>0
                            Jrg2=regiongrowing(Jse,x(k),y(k),th); 
                        end
                    end
                    % some morphological operations can be used here to
                    % correct the edges.
                    masks{it,i,is}=Jrg2;
        end
        %%
    end
end
end

