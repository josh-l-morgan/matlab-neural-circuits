function [param,finecorrparam]=pipeline_computefinelocalcorr(param)
%This function computes the actual peak of cross-correlation between local image regions
%Predicts and computes corresponding regions at the same time.
%For regions close to the edge, the algorithm attempts to move the region inside the image
%For use with the pipeline
%By Daniel Berger for MIT-BCS Seung, June 19 2009

%global tpatch;

lparam=param.computefinelocalcorr;

gridwidth=lparam.gridwidth;
gridheight=lparam.gridheight;
%patchwidth=lparam.patchwidth;
%patchheight=lparam.patchheight;
scaledpatchwidth=lparam.patchwidth/lparam.downscalefactor;
scaledpatchheight=lparam.patchheight/lparam.downscalefactor;
imagewidth=param.rawsize(2)/lparam.downscalefactor;
imageheight=param.rawsize(1)/lparam.downscalefactor;
midx=imagewidth/2;
midy=imageheight/2;

nrofslices=param.nrofslices;
nrofsrows=param.nrofrows;
nroftrows=param.nrofrows;
nrofscolumns=param.nrofcolumns;
nroftcolumns=param.nrofcolumns;
nroftslices=3; %previous, this, next

%corresponding coordinates of source-target pairs
%source is always a regular grid
% corresp_sx=zeros(nrofslices,nrofsrows,nrofscolumns,nroftslices,nroftrows,nroftcolumns,gridheight,gridwidth);
% corresp_sy=zeros(nrofslices,nrofsrows,nrofscolumns,nroftslices,nroftrows,nroftcolumns,gridheight,gridwidth);
%predcorresp_tx=zeros(nrofslices,nrofsrows,nrofscolumns,nroftslices,nroftrows,nroftcolumns,gridheight,gridwidth);
%predcorresp_ty=zeros(nrofslices,nrofsrows,nrofscolumns,nroftslices,nroftrows,nroftcolumns,gridheight,gridwidth);
%corresp_tphi=zeros(nrofslices,nrofsrows,nrofscolumns,nroftslices,nroftrows,nroftcolumns);

%Compute coordinates of source grid on source image, so that patches don't cross the image borders
sgridx=ones(gridheight,1)*[1:1:gridwidth]; %prepare grid as natural numbers
sgridy=[1:1:gridheight]'*ones(1,gridwidth); %sgridx';
sgridx=floor((sgridx-1)*scaledpatchwidth+scaledpatchwidth/2+1); %(sgridx-1)*((imagewidth-1)/(gridsize-1))+1; %scale grid to appropriate cordinates
sgridy=floor((sgridy-1)*scaledpatchheight+scaledpatchheight/2+1); %(sgridy-1)*((imageheight-1)/(gridsize-1))+1;

corresp_tx=zeros(gridheight,gridwidth);
corresp_ty=zeros(gridheight,gridwidth);

% corresp_sx=(param.predictcorr.corresp_sx-1)/lparam.downscalefactor+1;
% corresp_sy=(param.predictcorr.corresp_sy-1)/lparam.downscalefactor+1;
% corresp_tx=(param.predictcorr.corresp_tx-1)/lparam.downscalefactor+1;  <-- we'll compute this on the fly
% corresp_ty=(param.predictcorr.corresp_ty-1)/lparam.downscalefactor+1;
%corresp_tphi=param.predictcorr.corresp_tphi;

%The following variables will store the corresponding-point information
% corresp_csx=zeros(nrofslices,nrofsrows,nrofscolumns,nroftslices,nroftrows,nroftcolumns,gridheight,gridwidth);
% corresp_csy=zeros(nrofslices,nrofsrows,nrofscolumns,nroftslices,nroftrows,nroftcolumns,gridheight,gridwidth);
corresp_ctx=zeros(nrofslices,nrofsrows,nrofscolumns,nroftslices,nroftrows,nroftcolumns,gridheight,gridwidth);
corresp_cty=zeros(nrofslices,nrofsrows,nrofscolumns,nroftslices,nroftrows,nroftcolumns,gridheight,gridwidth);
corresp_cmax=zeros(nrofslices,nrofsrows,nrofscolumns,nroftslices,nroftrows,nroftcolumns,gridheight,gridwidth);
%corresp_cmean=zeros(nrofslices,nrofsrows,nrofscolumns,nroftslices,nroftrows,nroftcolumns,gridheight,gridwidth);
%corresp_cvar05=zeros(nrofslices,nrofsrows,nrofscolumns,nroftslices,nroftrows,nroftcolumns,gridheight,gridwidth);
corresp_cvar075=zeros(nrofslices,nrofsrows,nrofscolumns,nroftslices,nroftrows,nroftcolumns,gridheight,gridwidth);


%since the patch size is constant, we can pre-compute a patch grid which will 
%then be translated/rotated to be used for interp2 later
XI=ones(scaledpatchheight,1)*(1:scaledpatchwidth); XI=XI-scaledpatchwidth/2;
YI=(1:scaledpatchheight)'*ones(1,scaledpatchwidth); YI=YI-scaledpatchheight/2;

%A=[[0 1]; [1 0]]; A=repmat(A,32,32); %blank pattern

completetime=tic;

for ts=1:1:3
  %Set start and end slices according to target slice
  switch(ts)
    case 1 %previous slice
      beginslice=2;
      endslice=nrofslices;
    case 2 %same slice
      beginslice=1;
      endslice=nrofslices;
    case 3 %next slice
      beginslice=1;
      endslice=nrofslices-1;
  end;
  for slice=beginslice:1:endslice
  %for slice=168:1:endslice
    slicetimer=tic;
    for srow=1:1:nrofsrows
      for scol=1:1:nrofscolumns
        if (nrofsrows~=1)||(nroftrows~=1)||(nrofscolumns~=1)||(nroftcolumns~=1)||(ts~=2)  %simage not needed if nrofsrows,nrofscols,nroftrows,nroftcols=1 and ts=2.
          %Load source tile image
          filename=getrawfilename(param,slice,srow,scol);
          txt=sprintf('Loading source image %s ...',filename); disp(txt);
          simage=imread(filename);
          simage=double(simage)/255;
          if (lparam.downscalefactor~=1)
            simage=imresize(simage,1/lparam.downscalefactor);
          end;
        else
          filename=getrawfilename(param,slice,srow,scol);
          txt=sprintf('Skipping source image %s ...',filename); disp(txt);
        end;
        
        for trow=1:1:nroftrows
          for tcol=1:1:nroftcolumns
            
            if (ts==2)&&(trow==srow)&&(tcol==scol)
              %This is the exact same tile, so no computation is needed (tgrid=sgrid)
%               corresp_csx(slice,srow,scol,ts,trow,tcol,:,:)=corresp_sx(slice,srow,scol,ts,trow,tcol,:,:);
%               corresp_csy(slice,srow,scol,ts,trow,tcol,:,:)=corresp_sy(slice,srow,scol,ts,trow,tcol,:,:);
              corresp_ctx(slice,srow,scol,ts,trow,tcol,:,:)=sgridx; %corresp_sx(slice,srow,scol,ts,trow,tcol,:,:);
              corresp_cty(slice,srow,scol,ts,trow,tcol,:,:)=sgridy; %corresp_sy(slice,srow,scol,ts,trow,tcol,:,:);
              %We will leave cmax, cmean, cvar at zero for now.
            else
            
              %tstartimg=tic;
              %Load target tile image
              filename=getrawfilename(param,slice+ts-2,trow,tcol);
              txt=sprintf('Loading target image %s ...',filename); disp(txt);
              timage=imread(filename);
              timage=double(timage)/255;
              if (lparam.downscalefactor~=1)
                timage=imresize(timage,1/lparam.downscalefactor);
              end;
              
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              %predict source and target patch locations for this (slice,srow,scol,ts,trow,tcol)
              %source grid (corresp_sx/y) is always specified by sgridx and sgridy.
              %This computes corresp_tx, corresp_ty, corresp_tphi for the current loop only
              switch(ts)
                case 1 %with previous slice
                  smtx=squeeze(param.fititerativeaffine.affine(slice,srow,scol,:,:));
                  smtx(1,3)=smtx(1,3)/lparam.downscalefactor; %*param.downscale.scale;
                  smtx(2,3)=smtx(2,3)/lparam.downscalefactor; %*param.downscale.scale;
                  tmtx=squeeze(param.fititerativeaffine.affine(slice-1,trow,tcol,:,:));
                  tmtx(1,3)=tmtx(1,3)/lparam.downscalefactor; %*param.downscale.scale;
                  tmtx(2,3)=tmtx(2,3)/lparam.downscalefactor; %*param.downscale.scale;
                  rmtx=inv(tmtx)*smtx; %relative transformation matrix
                  
                  for gy=1:gridheight
                    parctx=zeros(1,gridwidth);
                    parcty=zeros(1,gridwidth);
                    for gx=1:gridwidth
                      tg=[sgridx(gy,gx)-midx sgridy(gy,gx)-midy 1]';
                      rtg=rmtx*tg;
                      parctx(gx)=rtg(1)+midx; %tgridx(gy,gx);
                      parcty(gx)=rtg(2)+midy; %tgridy(gy,gx);
                    end;
                    corresp_tx(gy,:)=parctx; %rtg(1)+midx; %tgridx(gy,gx);
                    corresp_ty(gy,:)=parcty; %rtg(2)+midy; %tgridy(gy,gx);
                  end;
                  corresp_tphi=atan2(rmtx(2,1),rmtx(1,1));
                  
                case 2 %with same slice
                  smtx=squeeze(param.fititerativeaffine.affine(slice,srow,scol,:,:));
                  smtx(1,3)=smtx(1,3)/lparam.downscalefactor;
                  smtx(2,3)=smtx(2,3)/lparam.downscalefactor;
                  tmtx=squeeze(param.fititerativeaffine.affine(slice,trow,tcol,:,:));
                  tmtx(1,3)=tmtx(1,3)/lparam.downscalefactor;
                  tmtx(2,3)=tmtx(2,3)/lparam.downscalefactor;
                  rmtx=inv(tmtx)*smtx; %relative transformation matrix
                  
                  for gy=1:1:gridheight
                    parctx=zeros(1,gridwidth);
                    parcty=zeros(1,gridwidth);
                    for gx=1:1:gridwidth
                      tg=[sgridx(gy,gx)-midx sgridy(gy,gx)-midy 1]';
                      rtg=rmtx*tg;
                      parctx(gx)=rtg(1)+midx;
                      parcty(gx)=rtg(2)+midy;
                    end;
                    corresp_tx(gy,:)=parctx;
                    corresp_ty(gy,:)=parcty;
                  end;
                  corresp_tphi=atan2(rmtx(2,1),rmtx(1,1));
                  
                case 3 %with next slice
                  smtx=squeeze(param.fititerativeaffine.affine(slice,srow,scol,:,:)); %We don't have relative transformations between arbitrary tiles
                  smtx(1,3)=smtx(1,3)/lparam.downscalefactor;
                  smtx(2,3)=smtx(2,3)/lparam.downscalefactor;
                  tmtx=squeeze(param.fititerativeaffine.affine(slice+1,trow,tcol,:,:));
                  tmtx(1,3)=tmtx(1,3)/lparam.downscalefactor;
                  tmtx(2,3)=tmtx(2,3)/lparam.downscalefactor;
                  rmtx=inv(tmtx)*smtx; %relative transformation matrix
                  
                  for gy=1:1:gridheight
                    parctx=zeros(1,gridwidth);
                    parcty=zeros(1,gridwidth);
                    for gx=1:1:gridwidth
                      tg=[sgridx(gy,gx)-midx sgridy(gy,gx)-midy 1]';
                      rtg=rmtx*tg;
                      parctx(gx)=rtg(1)+midx;
                      parcty(gx)=rtg(2)+midy;
                    end;
                    corresp_tx(gy,:)=parctx;
                    corresp_ty(gy,:)=parcty;
                  end;
                  corresp_tphi=atan2(rmtx(2,1),rmtx(1,1));
              end;
            
%             if (ts==2)&&(trow==srow)&&(tcol==scol)
%               %This is the exact same tile, so no computation is needed (tgrid=sgrid)
% %               corresp_csx(slice,srow,scol,ts,trow,tcol,:,:)=corresp_sx(slice,srow,scol,ts,trow,tcol,:,:);
% %               corresp_csy(slice,srow,scol,ts,trow,tcol,:,:)=corresp_sy(slice,srow,scol,ts,trow,tcol,:,:);
%               corresp_ctx(slice,srow,scol,ts,trow,tcol,:,:)=sgridx; %corresp_sx(slice,srow,scol,ts,trow,tcol,:,:);
%               corresp_cty(slice,srow,scol,ts,trow,tcol,:,:)=sgridy; %corresp_sy(slice,srow,scol,ts,trow,tcol,:,:);
%               %We will leave cmax, cmean, cvar at zero for now.
%             else
              for gy=1:1:gridheight
                for gx=1:1:gridwidth
                  %tstartpatch=tic;
                  %cut out image region from source image -
                  % this should always fit inside the image and rotation is 0
                  
%                   sx1=floor(corresp_sx(slice,srow,scol,ts,trow,tcol,gy,gx)-patchwidth/2);
%                   sx2=floor(sx1+patchwidth-1);
%                   sy1=floor(corresp_sy(slice,srow,scol,ts,trow,tcol,gy,gx)-patchheight/2);
%                   sy2=floor(sy1+patchheight-1);
%                   spatch=simage(sy1:sy2,sx1:sx2);
                  sx1=floor(sgridx(gy,gx)-scaledpatchwidth/2);
                  sx2=floor(sx1+scaledpatchwidth-1);
                  sy1=floor(sgridy(gy,gx)-scaledpatchheight/2);
                  sy2=floor(sy1+scaledpatchheight-1);
                  spatch=simage(sy1:sy2,sx1:sx2);
                  
                  %Cut out corresponding image region out of target image (rotated)
                  tgridx=corresp_tx(gy,gx);
                  tgridy=corresp_ty(gy,gx);
                  tphi=corresp_tphi;
                  
                  %Check whether the four patch corners are inside the image
                  cornersx=[-scaledpatchwidth/2 scaledpatchwidth/2 scaledpatchwidth/2 -scaledpatchwidth/2];
                  cornersy=[-scaledpatchheight/2 -scaledpatchheight/2 scaledpatchheight/2 scaledpatchheight/2];
                  rotcornersx=cos(tphi)*cornersx-sin(tphi)*cornersy;
                  rotcornersy=sin(tphi)*cornersx+cos(tphi)*cornersy;
                  rotcornersx=rotcornersx+tgridx;
                  rotcornersy=rotcornersy+tgridy;
                  
                  ctx=0; cty=0; %to store additional displacement
                  if (min(rotcornersx)>=1)&&(min(rotcornersy)>=1)&&(max(rotcornersx)<=size(timage,2))&&(max(rotcornersy)<=size(timage,1))
                    %This patch fits completely inside the target image.
                    tpatchfits=1;
                  else
                    %This patch does not fit inside the target image.
                    
                    %Check whether patch is completely outside of image
                    if (max(rotcornersx)<1)||(max(rotcornersy)<1)||(min(rotcornersx)>size(timage,2))||(min(rotcornersy)>size(timage,1))
                      %this patch is completely out of range.
                      tpatchfits=0;
                    else
                      %see whether the target patch can be fixed
                      %move target image part inwards; ctx and cty store that displacement
                      
                      if (min(rotcornersx)<1)
                        ctx=-min(rotcornersx)+1;
                      end;
                      if (min(rotcornersy)<1)
                        cty=-min(rotcornersy)+1;
                      end;
                      if (max(rotcornersx)>size(timage,2))
                        ctx=ctx+(size(timage,2)-max(rotcornersx))-1;
                      end;
                      if (max(rotcornersy)>size(timage,1))
                        cty=cty+(size(timage,1)-max(rotcornersy))-1;
                      end;
                      
                      %check whether there is still enough expected overlap
                      if (abs(ctx)<(scaledpatchwidth*(1-lparam.minoverlap)))&&(abs(cty)<(scaledpatchheight*(1-lparam.minoverlap)))
                        %overlap large enough, use this patch
                        tpatchfits=1;
                      else
                        %overlap too small
                        tpatchfits=0;
                      end;
                    end;
                  end;
                  
                  if tpatchfits==1
                    %use interp2 to read out rotated patch
                    XIr=cos(tphi)*XI-sin(tphi)*YI;
                    YIr=sin(tphi)*XI+cos(tphi)*YI;
                    XIr=XIr+tgridx+ctx;
                    YIr=YIr+tgridy+cty;

                    %compute rotated slice
%                     tic
%                     tpatch=interp2(timage,XIr,YIr,'linear'); %'nearest');
%                     toc

                    %CUT-OUT REGION FROM timage AND RELOCATE XIr, YIr
                    %tic
                    pminx=floor(min(min(XIr))); pmaxx=floor(max(max(XIr))+1);
                    pminy=floor(min(min(YIr))); pmaxy=floor(max(max(YIr))+1);
                    pimage=timage(pminy:pmaxy,pminx:pmaxx);
                    tpatch=interp2(pimage,XIr-pminx+1,YIr-pminy+1,'linear'); %'nearest');
                    %toc
                  end;
                  
                  %               if tpatchfits&&(isnan(tpatch(1,1)))
                  %                 disp('NAN!');
                  %               end;
                  
                  if lparam.showimages==1
                    if tpatchfits
                      figure(16);
                      imagesc(spatch);
                      axis equal;
                      colormap(gray);
                      title('Source patch');
                      
                      figure(17);
                      imagesc(tpatch);
                      axis equal;
                      colormap(gray);
                      if (ctx==0)&&(cty==0)
                        title('Target patch, FULL');
                      else
                        title('Target patch, PARTIAL');
                        hold on;
                        drawpoly([1 scaledpatchwidth scaledpatchwidth 1]-ctx,[1 1 scaledpatchheight scaledpatchheight]-cty,'r');
                        hold off;
                      end;
                      %                 else
                      %                   figure(17);
                      %                   imagesc(A);
                      %                   axis equal;
                      %                   colormap(gray);
                      %                   title('Target patch doesnt fit');
                    end;
                    pause(0.1);
                  end;
                  
                  if (tpatchfits==1)
                    %Do image filtering first if requested
                    if (lparam.dofiltering)
                      spatch=bandpass2(spatch,lparam.minppcyc,lparam.maxppcyc); %2000); for no highpass %=bandpass2(sprimage,40,80);
                      if lparam.showimages==1
                        figure(18);
                        imagesc(spatch);
                        axis equal;
                        colormap('gray');
                        title('Bandpass-filtered source patch');
                        axis square;
                      end;
                      
                      tpatch=bandpass2(tpatch,lparam.minppcyc,lparam.maxppcyc);
                      if lparam.showimages==1
                        figure(19);
                        imagesc(tpatch);
                        axis equal;
                        colormap('gray');
                        title('Bandpass-filtered target patch');
                        axis square;
                      end;
                    else %if no filtering: subtract mean to get rid of overlap size effect!
                      spatch=spatch-mean(mean(spatch));
                      if lparam.showimages==1
                        figure(18);
                        imagesc(spatch);
                        axis equal;
                        colormap('gray');
                        title('Bandpass-filtered source patch');
                        axis square;
                      end;
                      tpatch=tpatch-mean(mean(tpatch));
                      if lparam.showimages==1
                        figure(19);
                        imagesc(tpatch);
                        axis equal;
                        colormap('gray');
                        title('Bandpass-filtered target patch');
                        axis square;
                      end;                      
                    end;
                    
                    %Now compute the cross-correlation between spatch and tpatch
                    tic
                    if (ctx==0)&&(cty==0) %interior patch
                      transmap=convn_fast(spatch,flipdims(tpatch),'same');
                    else %Allow for larger displacements on boundary patches
                      transmap=convn_fast(spatch,flipdims(tpatch),'full');
                    end;
                    if lparam.showimages==1
                      figure(20);
                      imagesc(transmap);
                      pause(0.1);
                    end;
                    %toc
                    
                    %Compute maximum and other characteristics of cross-correlation result
                    middx=floor(size(transmap,2)/2+1);
                    middy=floor(size(transmap,1)/2+1);
                    [A,B]=max(transmap);
                    [C,D]=max(A);
                    transx=D-middx;  %X-translation of second image leading to maximal correlation
                    transy=B(D)-middy; %y-translation of second image leading to maximal correlation
                    transcorr=C;
                    %corresp_corr_tx(wafer,slice,srow,scol,ts,trow,tcol,gy,gx)=transx; %corresp_tx(wafer,slice,srow,scol,ts,trow,tcol,gy,gx)+transx;
                    %corresp_corr_ty(wafer,slice,srow,scol,ts,trow,tcol,gy,gx)=transy; %corresp_ty(wafer,slice,srow,scol,ts,trow,tcol,gy,gx)+transy;
                    if lparam.showimages==1
                      hold on;
                      %plot(middx+transx,middy+transy,'w*');
                      plot(middx+transx,middy+transy,'ko');
                      hold off;
                    end;
                    corresp_cmax(slice,srow,scol,ts,trow,tcol,gy,gx)=transcorr;
                    %corresp_cmean(slice,srow,scol,ts,trow,tcol,gy,gx)=mean(mean(transmap));
                    
                    %[x,y]=find(transmap>0.5*transcorr);
                    %xv=var(x); yv=var(y);
                    %corresp_cvar05(slice,srow,scol,ts,trow,tcol,gy,gx)=sqrt(xv*xv+yv*yv);
                    [x,y]=find(transmap>0.75*transcorr);
                    xv=var(x); yv=var(y);
                    corresp_cvar075(slice,srow,scol,ts,trow,tcol,gy,gx)=sqrt(xv*xv+yv*yv);
                    
                    if corresp_cvar075(slice,srow,scol,ts,trow,tcol,gy,gx)>100
                      col='r';
                    else
                      col='b';
                    end;
                    
                    %Compute center of overlapping region of source and target patch
                    smx1=max(1+transx,1); smx2=min(transx+scaledpatchwidth,scaledpatchwidth);
                    smy1=max(1+transy,1); smy2=min(transy+scaledpatchheight,scaledpatchheight);
                    smx=(smx1+smx2)/2;
                    smy=(smy1+smy2)/2;
                    tmx1=max(1-transx,1); tmx2=min(-transx+scaledpatchwidth,scaledpatchwidth);
                    tmy1=max(1-transy,1); tmy2=min(-transy+scaledpatchheight,scaledpatchheight);
                    tmx=(tmx1+tmx2)/2;
                    tmy=(tmy1+tmy2)/2;
                    rotmx=cos(tphi)*(tmx-scaledpatchwidth/2)-sin(tphi)*(tmy-scaledpatchheight/2)+tgridx;
                    rotmy=sin(tphi)*(tmx-scaledpatchwidth/2)+cos(tphi)*(tmy-scaledpatchheight/2)+tgridy;
                    
                    if lparam.showimages==1
                      figure(18);
                      hold on;
                      drawpoly([1 scaledpatchwidth scaledpatchwidth 1]+transx,[1 1 scaledpatchheight scaledpatchheight]+transy,col);
                      plot(smx,smy,'k*');
                      plot(smx,smy,'wo');
                      axis equal;
                      hold off;
                    end;
                    
                    if lparam.showimages==1
                      figure(19);
                      hold on;
                      drawpoly([1 scaledpatchwidth scaledpatchwidth 1]-transx,[1 1 scaledpatchheight scaledpatchheight]-transy,col);
                      plot(tmx,tmy,'k*');
                      plot(tmx,tmy,'wo');
                      axis equal;
                      hold off;
                    end;
                    
                    if lparam.showimages==1
                      figure(21);
                      imagesc(simage);
                      colormap('gray');
                      hold on;
                      drawpoly([sx1 sx2 sx2 sx1],[sy1 sy1 sy2 sy2],col);
                      plot(smx+sx1,smy+sy1,'k*');
                      plot(smx+sx1,smy+sy1,'wo');
                      hold off;
                      title('Source image part location');
                    end;
                    
                    if lparam.showimages==1
                      figure(22);
                      imagesc(timage);
                      colormap('gray');
                      hold on;
                      drawpoly(rotcornersx+ctx,rotcornersy+cty,col);
                      plot(rotmx+ctx,rotmy+cty,'k*');
                      plot(rotmx+ctx,rotmy+cty,'wo');
                      hold off;
                      title('Target image part location');
                    end;
                    
%                     corresp_csx(slice,srow,scol,ts,trow,tcol,gy,gx)=smx+sx1;
%                     corresp_csy(slice,srow,scol,ts,trow,tcol,gy,gx)=smy+sy1;
                    corresp_ctx(slice,srow,scol,ts,trow,tcol,gy,gx)=rotmx+ctx;
                    corresp_cty(slice,srow,scol,ts,trow,tcol,gy,gx)=rotmy+cty;
                    
                    %transcorr
                    %corresp_corr_vargt05(wafer,slice,srow,scol,ts,trow,tcol,gy,gx)
                    %corresp_corr_vargt075(wafer,slice,srow,scol,ts,trow,tcol,gy,gx)
                    
                    if (lparam.stopeach==1)
                      input('Press return for next');
                    end;
                    %if col=='r'
                    %  r=input('Is this correct?');
                    %end;
                  end;
                  %                 telapsedpatch=toc(tstartpatch);
                  %                 txt=sprintf('Patch duration: %f s',telapsedpatch);
                  %                 disp(txt);
                end;
              end;
            end;
%             telapsedimg=toc(tstartimg);
%             txt=sprintf('Image duration: %f s',telapsedimg);
%             disp(txt);
          end;
        end;
      end;
    end;
    lastslicetime=toc(slicetimer);
    txt=sprintf('Slice duration: %f s',lastslicetime);
    disp(txt);
  end;
end;

param.computefinelocalcorr.sgridx=sgridx;
param.computefinelocalcorr.sgridy=sgridy;
finecorrparam.computefinelocalcorr.sgridx=sgridx;
finecorrparam.computefinelocalcorr.sgridy=sgridy;
%param.computelocalcorr.corresp_csx=(corresp_csx-1)*lparam.downscalefactor+1; %the -1 +1 thing is needed so that scaling happens around 0
%param.computelocalcorr.corresp_csy=(corresp_csy-1)*lparam.downscalefactor+1;
finecorrparam.computefinelocalcorr.corresp_ctx=(corresp_ctx-1)*lparam.downscalefactor+1;
finecorrparam.computefinelocalcorr.corresp_cty=(corresp_cty-1)*lparam.downscalefactor+1;
finecorrparam.computefinelocalcorr.corresp_cmax=corresp_cmax;
%param.computelocalcorr.corresp_cmean=corresp_cmean;
%param.computelocalcorr.corresp_cvar05=corresp_cvar05; %*lparam.downscalefactor;
finecorrparam.computefinelocalcorr.corresp_cvar075=corresp_cvar075; %*lparam.downscalefactor;

lastcompletetime=toc(completetime);
hours=floor(lastcompletetime/3600);
mins=floor((lastcompletetime-hours*3600)/60);
secs=lastcompletetime-3600*hours-60*mins;
txt=sprintf('Complete duration: %f s (%d h, %d min, %f sec)',lastcompletetime,hours,mins,secs);
disp(txt);
