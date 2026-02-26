function pipeline_renderrigid(param)
%Renders out a stack of rigidly aligned images.
%For use with the pipeline
%By Daniel Berger fir MIT-BCS Seung, June 4th 2009

lparam=param.renderrigid;

imagewidth=param.scaledsize(2);
imageheight=param.scaledsize(1);

if (param.nrofrows==1)&&(param.nrofcolumns==1)
  if lparam.rendermode==1
    
    rigidmatrix=param.rigid.absrigidmatrix; %this has been estimated on either relative or absolute measures, but represented absolute
    
    %compute overall bounding box [minx,miny, maxx,maxy]
    firstslice=1;
    for slice=lparam.startslice:1:lparam.endslice
      %bounds=computeboundingbox(1,1,imagewidth,imageheight,squeeze(rigidmatrix(slice,:,:)));
      %x1=bounds(1); y1=bounds(2); x2=bounds(3); y2=bounds(4);
      bounds=computeboundingbox(-imagewidth/2,-imageheight/2,imagewidth/2,imageheight/2,squeeze(rigidmatrix(slice,1,1,:,:)));
      x1=bounds(1)+imagewidth/2; y1=bounds(2)+imageheight/2; x2=bounds(3)+imagewidth/2; y2=bounds(4)+imageheight/2;
      if firstslice
        firstslice=0;
        minx=x1; miny=y1; maxx=x2; maxy=y2;
      else
        minx=min([minx x1]); miny=min([miny y1]); maxx=max([maxx x2]); maxy=max([maxy y2]);
      end;
    end;
    minx=floor(minx); miny=floor(miny); maxx=floor(maxx+1); maxy=floor(maxy+1);
    
    %render out images
    twidth=maxx-minx+1; theight=maxy-miny+1;
    for slice=lparam.startslice:1:lparam.endslice
      name=sprintf(param.basescaledname,slice);
      filename=sprintf('%s%s',param.scaleddir,name);
      if exist(filename,'file')
        txt=sprintf('Loading image %s ...',filename);
        disp(txt);
        image=imread(filename);
        image=double(image)/255;
        
        image2=renderaffineimage(image,inv(squeeze(rigidmatrix(slice,1,1,:,:))),minx,miny,maxx,maxy,twidth,theight,'linear');
        name=sprintf(lparam.basename,slice);
        filename=sprintf('%s%s',lparam.targetdir,name);
        txt=sprintf(' Writing image %s ...',filename);
        disp(txt);
        imwrite(image2,filename,'png');
      end;
    end;
  else
    disp('ERROR: For single-tile stacks only rendermode 1 makes sense!');
  end;
else
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Multiple rows/columns
  
  rigidmatrix=param.rigid.absrigidmatrix; %this has been estimated on either relative or absolute measures, but represented absolute
  
  %compute overall bounding box [minx,miny, maxx,maxy]
  firstslice=1;
  for slice=lparam.startslice:1:lparam.endslice
    for row=1:1:param.nrofrows
      for column=1:1:param.nrofcolumns
        %bounds=computeboundingbox(1,1,imagewidth,imageheight,squeeze(rigidmatrix(slice,:,:)));
        %x1=bounds(1); y1=bounds(2); x2=bounds(3); y2=bounds(4);
        bounds=computeboundingbox(-imagewidth/2,-imageheight/2,imagewidth/2,imageheight/2,squeeze(rigidmatrix(slice,row,column,:,:)));
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
  
  %render out images
  twidth=maxx-minx+1; theight=maxy-miny+1;
  
  switch(lparam.rendermode)
    case 1   %rendermode 1: single tile per image
      disp('Using rendermode 1: single tile per image');
      for row=1:1:param.nrofrows
        for column=1:1:param.nrofcolumns
          for slice=lparam.startslice:1:lparam.endslice
            name=sprintf(param.basescaledname,slice,row,column);
            filename=sprintf('%s%s',param.scaleddir,name);
            if exist(filename,'file')
              txt=sprintf('Loading image %s ...',filename);
              disp(txt);
              image=imread(filename);
              image=double(image)/255;
              
              image2=renderaffineimage(image,inv(squeeze(rigidmatrix(slice,row,column,:,:))),minx,miny,maxx,maxy,twidth,theight,'linear');
              name=sprintf(lparam.basename,slice,row,column);
              filename=sprintf('%s%s',lparam.targetdir,name);
              txt=sprintf(' Writing image %s ...',filename);
              disp(txt);
              imwrite(image2,filename,'png');
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
            name=sprintf(param.basescaledname,slice,row,column);
            filename=sprintf('%s%s',param.scaleddir,name);
            if exist(filename,'file')
              txt=sprintf('Loading image %s ...',filename);
              disp(txt);
              image=imread(filename);
              image=double(image)/255;
              image2=renderaffineimage(image,inv(squeeze(rigidmatrix(slice,row,column,:,:))),minx,miny,maxx,maxy,twidth,theight,'linear');
              channels=[0 0 0]; channels(chan)=1;
              targetimage=copycoloredinto(targetimage,image2,channels);
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
      
      sineweightimage=generatesineweightimage(param.scaledsize(1),param.scaledsize(2));
      
      for slice=lparam.startslice:1:lparam.endslice
        targetweight=zeros(theight,twidth);
        targetimage=zeros(theight,twidth);
        for row=1:1:param.nrofrows
          for column=1:1:param.nrofcolumns
            name=sprintf(param.basescaledname,slice,row,column);
            filename=sprintf('%s%s',param.scaleddir,name);
            if exist(filename,'file')
              txt=sprintf('Loading image %s ...',filename);
              disp(txt);
              image=imread(filename);
              image=double(image)/255;
              
              image2=renderaffineimage(image,inv(squeeze(rigidmatrix(slice,row,column,:,:))),minx,miny,maxx,maxy,twidth,theight,'linear');
              imagew=renderaffineimage(sineweightimage,inv(squeeze(rigidmatrix(slice,row,column,:,:))),minx,miny,maxx,maxy,twidth,theight,'linear');
              targetimage=addimageinto(targetimage,image2.*imagew);
              targetweight=addimageinto(targetweight,imagew);              

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
