function pipeline_renderaffinetiled(param)
%Renders out a stack of affinely aligned images.
%Uses tiling to reduce memory usage during rendering
%For use with the pipeline
%By Daniel Berger fir MIT-BCS Seung, June 23rd 2009

lparam=param.renderaffine;

global sourceimage;
global source;
global sourcemask;

%prepare affine matrices
if lparam.usedownscaled==1
  imagewidth=param.scaledsize(2);
  imageheight=param.scaledsize(1);
  affinematrix=param.fititerativeaffine.affine; %this has been estimated on either relative or absolute measures, but represented absolute
  %Downscale affine matrix translations
  for slice=1:1:param.nrofslices
    for row=1:1:param.nrofrows
      for column=1:1:param.nrofcolumns
        affinematrix(slice,row,column,1,3)=affinematrix(slice,row,column,1,3)/param.downscale.scale;
        affinematrix(slice,row,column,2,3)=affinematrix(slice,row,column,2,3)/param.downscale.scale;
      end;
    end;
  end;
else
  imagewidth=param.rawsize(2);
  imageheight=param.rawsize(1);
  affinematrix=param.fititerativeaffine.affine; %this has been estimated on either relative or absolute measures, but represented absolute
end;


if (param.nrofrows==1)&&(param.nrofcolumns==1) %one-tile-per-slice-datasets
  if lparam.rendermode~=1
    disp('ERROR: For single-tile stacks only rendermode 1 makes sense!');
  else  
    %compute overall bounding box [minx,miny, maxx,maxy]
    firstslice=1;
%    for slice=lparam.startslice:1:lparam.endslice
    for slice=1:1:param.nrofslices %!! The overall bounding box (and target canvas size) should be computet over all slices!
      %this version is wrong.
      %bounds=computeboundingbox(1,1,imagewidth,imageheight,squeeze(rigidmatrix(slice,:,:)));
      %x1=bounds(1); y1=bounds(2); x2=bounds(3); y2=bounds(4);
      bounds=computeboundingbox(-imagewidth/2,-imageheight/2,imagewidth/2,imageheight/2,squeeze(affinematrix(slice,1,1,:,:)));
      x1=bounds(1)+imagewidth/2; y1=bounds(2)+imageheight/2; x2=bounds(3)+imagewidth/2; y2=bounds(4)+imageheight/2;
      %This is an alternative version in which the re-centering is done inside the function 
%       bounds=computeboundingbox_recentered(1,1,imagewidth,imageheight,squeeze(affinematrix(slice,1,1,:,:)));
%       x1=bounds(1); y1=bounds(2); x2=bounds(3); y2=bounds(4);
      
      if firstslice
        firstslice=0;
        minx=x1; miny=y1; maxx=x2; maxy=y2;
      else
        minx=min([minx x1]); miny=min([miny y1]); maxx=max([maxx x2]); maxy=max([maxy y2]);
      end;
    end;
    minx=floor(minx); miny=floor(miny); maxx=floor(maxx+1); maxy=floor(maxy+1);
    
    if lparam.crop==1
      minx=minx+lparam.crop_x1;
      miny=miny+lparam.crop_y1;
      maxx=minx+lparam.width-1;
      maxy=miny+lparam.height-1;
    end;
    
    %render out images
%     twidth=maxx-minx+1; theight=maxy-miny+1;
%     tilegridx1=1:lparam.tilesize:twidth;
%     tilegridy1=1:lparam.tilesize:theight;
    row=1; column=1;
    
    for slice=lparam.startslice:1:lparam.endslice
      if param.renderaffine.usedownscaled
        name=sprintf(param.basescaledname,slice);
        filename=sprintf('%s%s',param.scaleddir,name);
      else
        name=sprintf(param.baserawname,slice);
        filename=sprintf('%s%s',param.rawdir,name);  
      end;
      if exist(filename,'file')
        txt=sprintf('Loading image %s ...',filename);
        disp(txt);
        image=imread(filename);
        sourceimage=double(image)/255;
        
        tic;
        affinetransform_global_source_sourceimage_mask_minmax_tiled(inv(squeeze(affinematrix(slice,row,column,:,:))), minx,miny,maxx,maxy, lparam.tilesize, 'linear');
        toc
%         for tiley=tilegridy1
%           tileh=min(lparam.tilesize,theight-(tiley-1));
%           for tilex=tilegridx1
%             tilew=min(lparam.tilesize,twidth-(tilex-1));
% 
%             %image2=renderaffineimage(image,inv(squeeze(affinematrix(slice,1,1,:,:))),minx,miny,maxx,maxy,twidth,theight,'linear');
%           end;
%         end;
        name=sprintf(lparam.basename,slice);
        filename=sprintf('%s%s',lparam.targetdir,name);
        txt=sprintf(' Writing image %s ...',filename);
        disp(txt);
        imwrite(source,filename,'png');
      end;
    end;
  end;
else
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Multiple rows/columns
  
  %rigidmatrix=param.rigid.absrigidmatrix; %this has been estimated on either relative or absolute measures, but represented absolute
  
  %compute overall bounding box [minx,miny, maxx,maxy]
  firstslice=1;
  for slice=1:1:param.nrofslices
    for row=1:1:param.nrofrows
      for column=1:1:param.nrofcolumns
        %bounds=computeboundingbox(1,1,imagewidth,imageheight,squeeze(rigidmatrix(slice,:,:)));
        %x1=bounds(1); y1=bounds(2); x2=bounds(3); y2=bounds(4);
        bounds=computeboundingbox(-imagewidth/2,-imageheight/2,imagewidth/2,imageheight/2,squeeze(affinematrix(slice,row,column,:,:)));
        x1=bounds(1)+imagewidth/2; y1=bounds(2)+imageheight/2; x2=bounds(3)+imagewidth/2; y2=bounds(4)+imageheight/2;
        if firstslice
          firstslice=0;
          minx=x1; miny=y1; maxx=x2; maxy=y2;
        else
          minx=min([minx x1]); miny=min([miny y1]); maxx=max([maxx x2]); maxy=max([maxy y2]);
        end;
      end;
    end;
  end;
  minx=floor(minx); miny=floor(miny); maxx=floor(maxx+1); maxy=floor(maxy+1);
  
  if lparam.crop==1
    minx=minx+lparam.crop_x1;
    miny=miny+lparam.crop_y1;
    maxx=minx+lparam.crop_width-1;
    maxy=miny+lparam.crop_height-1;
  end;
  
  %render out images
  twidth=maxx-minx+1; theight=maxy-miny+1;
  
  switch(lparam.rendermode)
    case 1   %rendermode 1: single tile per image
      disp('Using rendermode 1: single tile per image');
      for row=1:1:param.nrofrows
        for column=1:1:param.nrofcolumns
          for slice=lparam.startslice:1:lparam.endslice
            if param.renderaffine.usedownscaled
              name=sprintf(param.basescaledname,slice,row,column);
              filename=sprintf('%s%s',param.scaleddir,name);
            else
              name=sprintf(param.baserawname,slice,row,column);
              filename=sprintf('%s%s',param.rawdir,name); 
            end;
            if exist(filename,'file')
              txt=sprintf('Loading image %s ...',filename);
              disp(txt);
              image=imread(filename);
              sourceimage=double(image)/255;
              
              %image2=renderaffineimage(image,inv(squeeze(affinematrix(slice,row,column,:,:))),minx,miny,maxx,maxy,twidth,theight,'linear');
              tic;
              affinetransform_global_source_sourceimage_mask_minmax_tiled(inv(squeeze(affinematrix(slice,row,column,:,:))),minx,miny,maxx,maxy,lparam.tilesize,'linear');
              toc
              name=sprintf(lparam.basename,slice,row,column);
              filename=sprintf('%s%s',lparam.targetdir,name);
              txt=sprintf(' Writing image %s ...',filename);
              disp(txt);
              imwrite(source,filename,'png');
            end;
          end;
        end;
      end;
    case 2   %rendermode 2: all tiles of one slice in one image, colored
      disp('Using rendermode 2: all tiles of one slice in one image, colored');
      for slice=lparam.startslice:1:lparam.endslice
        targetimage=zeros(theight,twidth,3);
        chan=1;
        for row=1:1:param.nrofrows
          for column=1:1:param.nrofcolumns
            if param.renderaffine.usedownscaled
              name=sprintf(param.basescaledname,slice,row,column);
              filename=sprintf('%s%s',param.scaleddir,name);
            else
              name=sprintf(param.baserawname,slice,row,column);
              filename=sprintf('%s%s',param.rawdir,name); 
            end;
            if exist(filename,'file')
              txt=sprintf('Loading image %s ...',filename);
              disp(txt);
              image=imread(filename);
              sourceimage=double(image)/255;
              
              %image2=renderaffineimage(image,inv(squeeze(affinematrix(slice,row,column,:,:))),minx,miny,maxx,maxy,twidth,theight,'linear');
              tic;
              affinetransform_global_source_sourceimage_mask_minmax_tiled(inv(squeeze(affinematrix(slice,row,column,:,:))),minx,miny,maxx,maxy,lparam.tilesize,'linear');
              toc
              channels=[0 0 0]; channels(chan)=1;
              targetimage=copycoloredinto(targetimage,source,channels);
              chan=chan+1; 
              if (chan>3)
                chan=1;
              end;
            end;
          end;
        end;
        name=sprintf(lparam.basename,slice,row,column);
        filename=sprintf('%s%s',lparam.targetdir,name);
        txt=sprintf(' Writing image %s ...',filename);
        disp(txt);
        imwrite(targetimage,filename,'png');
      end;
    case 3   %rendermode 3: all tiles of one slice in one image, blended
      disp('Using rendermode 3: all tiles of one slice in one image, blended');
      if param.renderaffine.usedownscaled
        sineweightimage=generatesineweightimage(param.scaledsize(1),param.scaledsize(2));
      else
        sineweightimage=generatesineweightimage(param.rawsize(1),param.rawsize(2));
      end;
      
      for slice=lparam.startslice:1:lparam.endslice
        targetweight=zeros(theight,twidth);
        targetimage=zeros(theight,twidth);
        for row=1:1:param.nrofrows
          for column=1:1:param.nrofcolumns
            if param.renderaffine.usedownscaled
              name=sprintf(param.basescaledname,slice,row,column);
              filename=sprintf('%s%s',param.scaleddir,name);
            else
              name=sprintf(param.baserawname,slice,row,column);
              filename=sprintf('%s%s',param.rawdir,name); 
            end;
            if exist(filename,'file')
              txt=sprintf('Loading image %s ...',filename);
              disp(txt);
              image=imread(filename);
              sourceimage=double(image)/255;
              clear image;
              
              %image2=renderaffineimage(image,inv(squeeze(affinematrix(slice,row,column,:,:))),minx,miny,maxx,maxy,twidth,theight,'linear');
              %imagew=renderaffineimage(sineweightimage,inv(squeeze(affinematrix(slice,row,column,:,:))),minx,miny,maxx,maxy,twidth,theight,'linear');
              affinetransform_global_source_sourceimage_mask_minmax_tiled(inv(squeeze(affinematrix(slice,row,column,:,:))),minx,miny,maxx,maxy,lparam.tilesize,'linear');
              image2=source; sourceimage=sineweightimage;
              affinetransform_global_source_sourceimage_mask_minmax_tiled(inv(squeeze(affinematrix(slice,row,column,:,:))),minx,miny,maxx,maxy,lparam.tilesize,'linear');
              targetimage=addimageinto(targetimage,image2.*source);
              targetweight=addimageinto(targetweight,source);              
            end;
          end;
        end;
        
        %normalize image
        for th=1:1:theight
          for tw=1:1:twidth
            if targetweight(th,tw)>0
              targetimage(th,tw)=targetimage(th,tw)/targetweight(th,tw);
            end;
          end;
        end;
        name=sprintf(lparam.basename,slice,row,column);
        filename=sprintf('%s%s',lparam.targetdir,name);
        txt=sprintf(' Writing image %s ...',filename);
        disp(txt);
        imwrite(targetimage,filename,'png');
      end;
    otherwise
      disp('ERROR: Unknown rendermode.');
  end;
end;
