function affinetransform_global_source_sourceimage_mask_minmax_tiled(A, xmin,ymin,xmax,ymax, tilesize, method)
%Computes an affinely transformed image from image sourceimage and stores it in source.
%uses an evaluation grid specified by [xmin,ymin]-[xmax,ymax] (which should be integers).
%method is the method for interp2, for example 'linear'.
%This version of the version does interp2 in tiles, to reduce memory consumption.
%by Daniel Berger for MIT-BCS Seung, March 21 2009 (updated May 5 2009)

global sourceimage; %this is the source image
global source; %this is the target image (called source for further processing)
global sourcemask;

border=10;

%midx=size(sourceimage,2)/2; midy=size(sourceimage,1)/2;
%actually, what we need here, is the center of the sourceimage in  xmin..xmax, ymin..ymax coordinates
midx=size(sourceimage,2)/2; midy=size(sourceimage,1)/2;

targetxsize=xmax-xmin+1;
targetysize=ymax-ymin+1;

source=zeros(targetysize,targetxsize);
sourcemask=zeros(targetysize,targetxsize);

%compute bounding box for sourceimage
% Q=inv(A)*[[1 size(sourceimage,2) size(sourceimage,2) 1]; [1 1 size(sourceimage,1) size(sourceimage,1)]; [1 1 1 1]];
% bbminx=min(Q(1,:));
% bbmaxx=max(Q(1,:));
% bbminy=min(Q(2,:));
% bbmaxy=max(Q(2,:));
  bounds=computeboundingbox_recentered(1,1,size(sourceimage,2),size(sourceimage,1),inv(A));
  bbminx=bounds(1); bbminy=bounds(2); bbmaxx=bounds(3); bbmaxy=bounds(4);
  %bbminx=bounds(1)-xmin; bbminy=bounds(2)-ymin; bbmaxx=bounds(3)-xmin; bbmaxy=bounds(4)-ymin; %this would be for a bounding box in target image space

for y=ymin:tilesize:ymax
  y
  for x=xmin:tilesize:xmax
    %x
    tilexs=min(tilesize,xmax-x+1);
    tileys=min(tilesize,ymax-y+1);
    
    if (boundingboxoverlap(bbminx,bbminy,bbmaxx,bbmaxy,x,y,x+tilexs-1,y+tileys-1))
      
      %generate tile's coordinate grid
      XI=ones(tileys,1)*(x:x+tilexs-1);
      YI=(y:y+tileys-1)'*ones(1,tilexs);
      
      %transform coordinate grid
      XI=XI-midx; YI=YI-midy;
      XIa=XI; YIa=YI; %initialize to same size for increased speed
      for ly=1:1:size(XI,1)
        for lx=1:1:size(XI,2)
          XIa(ly,lx)=XI(ly,lx)*A(1,1)+YI(ly,lx)*A(1,2)+A(1,3);
          YIa(ly,lx)=XI(ly,lx)*A(2,1)+YI(ly,lx)*A(2,2)+A(2,3);
        end;
      end;
      XIa=XIa+midx; YIa=YIa+midy;
      
      %transform coordinate grid
%       XIa=XI*A(1,1)+YI*A(1,2)+A(1,3);
%       YIa=XI*A(2,1)+YI*A(2,2)+A(2,3);
      
      %find region of interest on sourceimage, clip, crop, adapt
      %transformation matrices
      ixmin=floor(min(min(XIa)))-border; ixmax=floor(max(max(XIa)+1))+border;
      iymin=floor(min(min(YIa)))-border; iymax=floor(max(max(YIa)+1))+border;
      if (ixmax>0)&&(iymax>0)&&(ixmin<size(sourceimage,2))&&(iymin<size(sourceimage,2))
        cixmin=max(ixmin,1); cixmax=min(ixmax,size(sourceimage,2));
        ciymin=max(iymin,1); ciymax=min(iymax,size(sourceimage,1));
        
        if ((cixmax-cixmin)>2)&&((ciymax-ciymin)>2)
          
          sourceimagepatch=sourceimage(ciymin:ciymax,cixmin:cixmax);
%           figure(60);
%           imagesc(sourceimagepatch);
%           colormap(gray);
%           title('Sourceimage patch');
          
          XIa=XIa+(1-cixmin);
          YIa=YIa+(1-ciymin);
          
          %compute affinely transformed tiles of image and mask
          sourcetile=interp2(sourceimagepatch,XIa,YIa,method);
          sourcemasktile=interp2(sourceimagepatch*0,XIa,YIa,method);
          
%           figure(61);
%           imagesc(sourcetile);
%           colormap(gray);
%           title('Computed tile');
%           input('Press Return ...');
        else
          sourcetile=ones(tileys,tilexs)*NaN;
          sourcemasktile=sourcetile;
        end;
      else
        sourcetile=ones(tileys,tilexs)*NaN;
        sourcemasktile=sourcetile;  
      end;
    
    else
      sourcetile=ones(tileys,tilexs)*NaN;
      sourcemasktile=sourcetile;
    end;
    
    %copy tiles into target images (source and sourcemask)
    ty=y-ymin+1; tx=x-xmin+1;
    source(ty:ty+tileys-1,tx:tx+tilexs-1)=sourcetile;
    sourcemask(ty:ty+tileys-1,tx:tx+tilexs-1)=sourcemasktile;
  end;
end;