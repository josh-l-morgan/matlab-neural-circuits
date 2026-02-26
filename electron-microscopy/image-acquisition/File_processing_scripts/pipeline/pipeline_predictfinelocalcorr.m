function [param,finecorrparam]=pipeline_predictfinelocalcorr(param)
%A function which computes where corresponding points in neighboring slices/tiles
%should be, based on rigid alignment
%For use with the pipeline
%By Daniel Berger for MIT-BCS Seung, June 6th 2009

lparam=param.predictfinelocalcorr;

gridwidth=lparam.gridwidth;
gridheight=lparam.gridheight;
patchwidth=lparam.patchwidth;
patchheight=lparam.patchheight;
imagewidth=param.rawsize(2);
imageheight=param.rawsize(1);
midx=imagewidth/2;
midy=imageheight/2;

nrofslices=param.nrofslices;
nrofsrows=param.nrofrows;
nroftrows=param.nrofrows;
nrofscolumns=param.nrofcolumns;
nroftcolumns=param.nrofcolumns;
nroftslices=3; %previous, this, next
param.predictcorr.nroftslices=nroftslices;
%meanrelrotang=squeeze(mean(mean(param.rigid.relrot,3),2));

%corresponding coordinates of source-target pairs
%source is always a regular grid
% corresp_sx=zeros(nrofslices,nrofsrows,nrofscolumns,nroftslices,nroftrows,nroftcolumns,gridheight,gridwidth);
% corresp_sy=zeros(nrofslices,nrofsrows,nrofscolumns,nroftslices,nroftrows,nroftcolumns,gridheight,gridwidth);
corresp_tx=zeros(nrofslices,nrofsrows,nrofscolumns,nroftslices,nroftrows,nroftcolumns,gridheight,gridwidth);
corresp_ty=zeros(nrofslices,nrofsrows,nrofscolumns,nroftslices,nroftrows,nroftcolumns,gridheight,gridwidth);
corresp_tphi=zeros(nrofslices,nrofsrows,nrofscolumns,nroftslices,nroftrows,nroftcolumns);

%Compute coordinates of source grid on source image, so that patches don't cross the image borders
sgridx=ones(gridheight,1)*[1:1:gridwidth]; %prepare grid as natural numbers
sgridy=[1:1:gridheight]'*ones(1,gridwidth); %sgridx';
sgridx=floor((sgridx-1)*patchwidth+patchwidth/2+1); %(sgridx-1)*((imagewidth-1)/(gridsize-1))+1; %scale grid to appropriate cordinates
sgridy=floor((sgridy-1)*patchheight+patchheight/2+1); %(sgridy-1)*((imageheight-1)/(gridsize-1))+1;

% disp('Initializing source coordinates...');
% %initialize source coordinate matrices
% for slice=1:1:nrofslices
%   for srow=1:1:nrofsrows
%     for scol=1:1:nrofscolumns
%       for ts=1:1:nroftslices
%         for trow=1:1:nroftrows
%           for tcol=1:1:nroftcolumns
%             for gy=1:1:gridheight
%               for gx=1:1:gridwidth
%                 corresp_sx(slice,srow,scol,ts,trow,tcol,gy,gx)=sgridx(gy,gx);
%                 corresp_sy(slice,srow,scol,ts,trow,tcol,gy,gx)=sgridy(gy,gx);
%               end;
%             end;
%           end;
%         end;
%       end;
%     end;
%   end;
% end;


%compute previous-slice target coordinate matrix
%this is now based on the absolute rigid transformation matrices for all tiles
%ATTENTION: this computes the grid coordinates in the LOCAL coordinate frame of the PREVIOUS-slice tile
%(where the corresponding points of the grid defined axis-parallel in the given source tile are in the previous slice)

%matlabpool open



ts=1; %from slice to previous slice
disp('Computing corresponding points to previous slice...');
for slice=2:1:nrofslices
  tic
  txt=sprintf('-> Prev, Slice %i/%i',slice,nrofslices);
  disp(txt);
  for srow=1:1:param.nrofrows
    for scol=1:1:param.nrofcolumns
      for trow=1:1:param.nrofrows
        for tcol=1:1:param.nrofcolumns
          
          smtx=squeeze(param.fititerativeaffine.affine(slice,srow,scol,:,:));
          %smtx=squeeze(param.rigid.absrigidmatrix(slice,srow,scol,:,:)); %use absolute because we don't have relative transformations between arbitrary tiles
          smtx(1,3)=smtx(1,3)*param.downscale.scale;
          smtx(2,3)=smtx(2,3)*param.downscale.scale;
          tmtx=squeeze(param.fititerativeaffine.affine(slice-1,trow,tcol,:,:));
          %tmtx=squeeze(param.rigid.absrigidmatrix(slice-1,trow,tcol,:,:));
          tmtx(1,3)=tmtx(1,3)*param.downscale.scale;
          tmtx(2,3)=tmtx(2,3)*param.downscale.scale;
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
            corresp_tx(slice,srow,scol,ts,trow,tcol,gy,:)=parctx; %rtg(1)+midx; %tgridx(gy,gx);
            corresp_ty(slice,srow,scol,ts,trow,tcol,gy,:)=parcty; %rtg(2)+midy; %tgridy(gy,gx);
          end;

          corresp_tphi(slice,srow,scol,ts,trow,tcol)=atan2(rmtx(2,1),rmtx(1,1)); %atan2(sin(1),cos(1));
        end;
      end;
    end;
  end;
  toc
end;


%compute same-slice target coordinate matrix
ts=2; %from slice to next slice
disp('Computing corresponding points between tiles in the same slice...');
for slice=1:1:nrofslices
  txt=sprintf('-> Same, Slice %i/%i',slice,nrofslices);
  disp(txt);
  for srow=1:1:param.nrofrows
    for scol=1:1:param.nrofcolumns
      for trow=1:1:param.nrofrows
        for tcol=1:1:param.nrofcolumns
          smtx=squeeze(param.fititerativeaffine.affine(slice,srow,scol,:,:));
          %smtx=squeeze(param.rigid.absrigidmatrix(slice,srow,scol,:,:));
          smtx(1,3)=smtx(1,3)*param.downscale.scale;
          smtx(2,3)=smtx(2,3)*param.downscale.scale;
          tmtx=squeeze(param.fititerativeaffine.affine(slice,trow,tcol,:,:));
          %tmtx=squeeze(param.rigid.absrigidmatrix(slice,trow,tcol,:,:));
          tmtx(1,3)=tmtx(1,3)*param.downscale.scale;
          tmtx(2,3)=tmtx(2,3)*param.downscale.scale;
          rmtx=inv(tmtx)*smtx; %relative transformation matrix
          
          for gy=1:1:gridheight
            parctx=zeros(1,gridwidth);
            parcty=zeros(1,gridwidth);
            for gx=1:1:gridwidth
              tg=[sgridx(gy,gx)-midx sgridy(gy,gx)-midy 1]';
              rtg=rmtx*tg;
              parctx(gx)=rtg(1)+midx; %tgridx(gy,gx);
              parcty(gx)=rtg(2)+midy; %tgridy(gy,gx);
            end;
            corresp_tx(slice,srow,scol,ts,trow,tcol,gy,:)=parctx; %tgridx(gy,gx);
            corresp_ty(slice,srow,scol,ts,trow,tcol,gy,:)=parcty; %tgridy(gy,gx);
          end;
          corresp_tphi(slice,srow,scol,ts,trow,tcol)=atan2(rmtx(2,1),rmtx(1,1)); %atan2(sin(1),cos(1));
        end;
      end;
    end;
  end;
end;


%compute next-slice target coordinate matrix
ts=3; %from slice to next slice
disp('Computing corresponding points to next slice...');
for slice=1:1:nrofslices-1
  txt=sprintf('-> Next, Slice %i/%i',slice,nrofslices);
  disp(txt);
  for srow=1:1:param.nrofrows
    for scol=1:1:param.nrofcolumns
      for trow=1:1:param.nrofrows
        for tcol=1:1:param.nrofcolumns
          smtx=squeeze(param.fititerativeaffine.affine(slice,srow,scol,:,:));
          %smtx=squeeze(param.rigid.absrigidmatrix(slice,srow,scol,:,:)); %We don't have relative transformations between arbitrary tiles
          smtx(1,3)=smtx(1,3)*param.downscale.scale;
          smtx(2,3)=smtx(2,3)*param.downscale.scale;
          tmtx=squeeze(param.fititerativeaffine.affine(slice+1,trow,tcol,:,:));
          %tmtx=squeeze(param.rigid.absrigidmatrix(slice+1,trow,tcol,:,:));
          tmtx(1,3)=tmtx(1,3)*param.downscale.scale;
          tmtx(2,3)=tmtx(2,3)*param.downscale.scale;
          rmtx=inv(tmtx)*smtx; %relative transformation matrix
          
          for gy=1:1:gridheight
            parctx=zeros(1,gridwidth);
            parcty=zeros(1,gridwidth);
            for gx=1:1:gridwidth
              tg=[sgridx(gy,gx)-midx sgridy(gy,gx)-midy 1]';
              rtg=rmtx*tg;
              parctx(gx)=rtg(1)+midx; %tgridx(gy,gx);
              parcty(gx)=rtg(2)+midy; %tgridy(gy,gx);
            end;
            corresp_tx(slice,srow,scol,ts,trow,tcol,gy,:)=parctx; %tgridx(gy,gx);
            corresp_ty(slice,srow,scol,ts,trow,tcol,gy,:)=parcty; %tgridy(gy,gx);
          end;
          corresp_tphi(slice,srow,scol,ts,trow,tcol)=atan2(rmtx(2,1),rmtx(1,1)); %atan2(sin(1),cos(1));
        end;
      end;
    end;
  end;
end;

%storing the source coordinates for each slice, srow, scol, ts, trow, tcol is a waste of memory, esp. for so large grids
finecorrparam.sgridx=sgridx;
finecorrparam.sgridy=sgridy;
% param.predictcorr.corresp_sx=corresp_sx;
% param.predictcorr.corresp_sy=corresp_sy;
finecorrparam.corresp_tx=corresp_tx;
finecorrparam.corresp_ty=corresp_ty;
finecorrparam.corresp_tphi=corresp_tphi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%% Plot for validation this->prev
% slice=168; %28
% ts=1;  %target slice: 1:previous, 2:same, 3:next
% 
% srow=1; scol=1; trow=1; tcol=1;
% %sradang=absolute_angle(wafer,slice,srow,scol)*pi/180;
% %tradang=absolute_angle(prevslicewafer,prevsliceslice,trow,tcol)*pi/180;
% sgridx=squeeze(corresp_sx(slice,srow,scol,ts,trow,tcol,:,:));
% sgridy=squeeze(corresp_sy(slice,srow,scol,ts,trow,tcol,:,:));
% tgridx=squeeze(corresp_tx(slice,srow,scol,ts,trow,tcol,:,:));
% tgridy=squeeze(corresp_ty(slice,srow,scol,ts,trow,tcol,:,:));
% name=sprintf(param.basescaledname,slice,srow,scol); filename=sprintf('%s%s',param.scaleddir,name);
% drawimgplusgrid2(30,filename,sgridx/param.downscale.scale,sgridy/param.downscale.scale,0,patchwidth/param.downscale.scale,patchheight/param.downscale.scale);
% txt=sprintf('Source: Slice %d, Row %d, Column %d',slice,srow,scol);
% title(txt);
% name=sprintf(param.basescaledname,slice-1,trow,tcol); filename=sprintf('%s%s',param.scaleddir,name);
% drawimgplusgrid2(31,filename,tgridx/param.downscale.scale,tgridy/param.downscale.scale,corresp_tphi(slice,srow,scol,ts,trow,tcol),patchwidth/param.downscale.scale,patchheight/param.downscale.scale);
% txt=sprintf('Target: Slice %d, Row %d, Column %d',slice-1,trow,tcol);
% title(txt);
% 
% trow=1; tcol=2;
% tgridx=squeeze(corresp_tx(slice,srow,scol,ts,trow,tcol,:,:));
% tgridy=squeeze(corresp_ty(slice,srow,scol,ts,trow,tcol,:,:));
% name=sprintf(param.basescaledname,slice-1,trow,tcol); filename=sprintf('%s%s',param.scaleddir,name);
% drawimgplusgrid2(32,filename,tgridx/param.downscale.scale,tgridy/param.downscale.scale,corresp_tphi(slice,srow,scol,ts,trow,tcol),patchwidth/param.downscale.scale,patchheight/param.downscale.scale);
% txt=sprintf('Target: Slice %d, Row %d, Column %d',slice-1,trow,tcol);
% title(txt);
% 
% trow=2; tcol=1;
% tgridx=squeeze(corresp_tx(slice,srow,scol,ts,trow,tcol,:,:));
% tgridy=squeeze(corresp_ty(slice,srow,scol,ts,trow,tcol,:,:));
% name=sprintf(param.basescaledname,slice-1,trow,tcol); filename=sprintf('%s%s',param.scaleddir,name);
% drawimgplusgrid2(33,filename,tgridx/param.downscale.scale,tgridy/param.downscale.scale,corresp_tphi(slice,srow,scol,ts,trow,tcol),patchwidth/param.downscale.scale,patchheight/param.downscale.scale);
% txt=sprintf('Target: Slice %d, Row %d, Column %d',slice-1,trow,tcol);
% title(txt);
% 
% trow=2; tcol=2;
% tgridx=squeeze(corresp_tx(slice,srow,scol,ts,trow,tcol,:,:));
% tgridy=squeeze(corresp_ty(slice,srow,scol,ts,trow,tcol,:,:));
% name=sprintf(param.basescaledname,slice-1,trow,tcol); filename=sprintf('%s%s',param.scaleddir,name);
% drawimgplusgrid2(34,filename,tgridx/param.downscale.scale,tgridy/param.downscale.scale,corresp_tphi(slice,srow,scol,ts,trow,tcol),patchwidth/param.downscale.scale,patchheight/param.downscale.scale);
% txt=sprintf('Target: Slice %d, Row %d, Column %d',slice-1,trow,tcol);
% title(txt);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%% Plot for validation this->this
% slice=2; %28
% ts=2;  %target slice: 1:previous, 2:same, 3:next
% 
% srow=1; scol=1; trow=1; tcol=1;
% %sradang=absolute_angle(wafer,slice,srow,scol)*pi/180;
% %tradang=absolute_angle(prevslicewafer,prevsliceslice,trow,tcol)*pi/180;
% sgridx=squeeze(corresp_sx(slice,srow,scol,ts,trow,tcol,:,:));
% sgridy=squeeze(corresp_sy(slice,srow,scol,ts,trow,tcol,:,:));
% tgridx=squeeze(corresp_tx(slice,srow,scol,ts,trow,tcol,:,:));
% tgridy=squeeze(corresp_ty(slice,srow,scol,ts,trow,tcol,:,:));
% name=sprintf(param.basescaledname,slice,srow,scol); filename=sprintf('%s%s',param.scaleddir,name);
% drawimgplusgrid2(40,filename,sgridx/param.downscale.scale,sgridy/param.downscale.scale,0,patchwidth/param.downscale.scale,patchheight/param.downscale.scale);
% txt=sprintf('Source: Slice %d, Row %d, Column %d',slice,srow,scol);
% title(txt);
% name=sprintf(param.basescaledname,slice,trow,tcol); filename=sprintf('%s%s',param.scaleddir,name);
% drawimgplusgrid2(41,filename,tgridx/param.downscale.scale,tgridy/param.downscale.scale,corresp_tphi(slice,srow,scol,ts,trow,tcol),patchwidth/param.downscale.scale,patchheight/param.downscale.scale);
% txt=sprintf('Target: Slice %d, Row %d, Column %d',slice,trow,tcol);
% title(txt);
% 
% trow=1; tcol=2;
% tgridx=squeeze(corresp_tx(slice,srow,scol,ts,trow,tcol,:,:));
% tgridy=squeeze(corresp_ty(slice,srow,scol,ts,trow,tcol,:,:));
% name=sprintf(param.basescaledname,slice,trow,tcol); filename=sprintf('%s%s',param.scaleddir,name);
% drawimgplusgrid2(42,filename,tgridx/param.downscale.scale,tgridy/param.downscale.scale,corresp_tphi(slice,srow,scol,ts,trow,tcol),patchwidth/param.downscale.scale,patchheight/param.downscale.scale);
% txt=sprintf('Target: Slice %d, Row %d, Column %d',slice,trow,tcol);
% title(txt);
% 
% trow=2; tcol=1;
% tgridx=squeeze(corresp_tx(slice,srow,scol,ts,trow,tcol,:,:));
% tgridy=squeeze(corresp_ty(slice,srow,scol,ts,trow,tcol,:,:));
% name=sprintf(param.basescaledname,slice,trow,tcol); filename=sprintf('%s%s',param.scaleddir,name);
% drawimgplusgrid2(43,filename,tgridx/param.downscale.scale,tgridy/param.downscale.scale,corresp_tphi(slice,srow,scol,ts,trow,tcol),patchwidth/param.downscale.scale,patchheight/param.downscale.scale);
% txt=sprintf('Target: Slice %d, Row %d, Column %d',slice,trow,tcol);
% title(txt);
% 
% trow=2; tcol=2;
% tgridx=squeeze(corresp_tx(slice,srow,scol,ts,trow,tcol,:,:));
% tgridy=squeeze(corresp_ty(slice,srow,scol,ts,trow,tcol,:,:));
% name=sprintf(param.basescaledname,slice,trow,tcol); filename=sprintf('%s%s',param.scaleddir,name);
% drawimgplusgrid2(44,filename,tgridx/param.downscale.scale,tgridy/param.downscale.scale,corresp_tphi(slice,srow,scol,ts,trow,tcol),patchwidth/param.downscale.scale,patchheight/param.downscale.scale);
% txt=sprintf('Target: Slice %d, Row %d, Column %d',slice,trow,tcol);
% title(txt);
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%% Plot for validation this->next
% slice=167; %28
% ts=3;  %target slice: 1:previous, 2:same, 3:next
% 
% srow=1; scol=1; trow=1; tcol=1;
% %sradang=absolute_angle(wafer,slice,srow,scol)*pi/180;
% %tradang=absolute_angle(prevslicewafer,prevsliceslice,trow,tcol)*pi/180;
% sgridx=squeeze(corresp_sx(slice,srow,scol,ts,trow,tcol,:,:));
% sgridy=squeeze(corresp_sy(slice,srow,scol,ts,trow,tcol,:,:));
% tgridx=squeeze(corresp_tx(slice,srow,scol,ts,trow,tcol,:,:));
% tgridy=squeeze(corresp_ty(slice,srow,scol,ts,trow,tcol,:,:));
% name=sprintf(param.basescaledname,slice,srow,scol); filename=sprintf('%s%s',param.scaleddir,name);
% drawimgplusgrid2(50,filename,sgridx/param.downscale.scale,sgridy/param.downscale.scale,0,patchwidth/param.downscale.scale,patchheight/param.downscale.scale);
% txt=sprintf('Source: Slice %d, Row %d, Column %d',slice,srow,scol);
% title(txt);
% name=sprintf(param.basescaledname,slice+1,trow,tcol); filename=sprintf('%s%s',param.scaleddir,name);
% drawimgplusgrid2(51,filename,tgridx/param.downscale.scale,tgridy/param.downscale.scale,corresp_tphi(slice,srow,scol,ts,trow,tcol),patchwidth/param.downscale.scale,patchheight/param.downscale.scale);
% txt=sprintf('Target: Slice %d, Row %d, Column %d',slice+1,trow,tcol);
% title(txt);
% 
% trow=1; tcol=2;
% tgridx=squeeze(corresp_tx(slice,srow,scol,ts,trow,tcol,:,:));
% tgridy=squeeze(corresp_ty(slice,srow,scol,ts,trow,tcol,:,:));
% name=sprintf(param.basescaledname,slice+1,trow,tcol); filename=sprintf('%s%s',param.scaleddir,name);
% drawimgplusgrid2(52,filename,tgridx/param.downscale.scale,tgridy/param.downscale.scale,corresp_tphi(slice,srow,scol,ts,trow,tcol),patchwidth/param.downscale.scale,patchheight/param.downscale.scale);
% txt=sprintf('Target: Slice %d, Row %d, Column %d',slice+1,trow,tcol);
% title(txt);
% 
% trow=2; tcol=1;
% tgridx=squeeze(corresp_tx(slice,srow,scol,ts,trow,tcol,:,:));
% tgridy=squeeze(corresp_ty(slice,srow,scol,ts,trow,tcol,:,:));
% name=sprintf(param.basescaledname,slice+1,trow,tcol); filename=sprintf('%s%s',param.scaleddir,name);
% drawimgplusgrid2(53,filename,tgridx/param.downscale.scale,tgridy/param.downscale.scale,corresp_tphi(slice,srow,scol,ts,trow,tcol),patchwidth/param.downscale.scale,patchheight/param.downscale.scale);
% txt=sprintf('Target: Slice %d, Row %d, Column %d',slice+1,trow,tcol);
% title(txt);
% 
% trow=2; tcol=2;
% tgridx=squeeze(corresp_tx(slice,srow,scol,ts,trow,tcol,:,:));
% tgridy=squeeze(corresp_ty(slice,srow,scol,ts,trow,tcol,:,:));
% name=sprintf(param.basescaledname,slice+1,trow,tcol); filename=sprintf('%s%s',param.scaleddir,name);
% drawimgplusgrid2(54,filename,tgridx/param.downscale.scale,tgridy/param.downscale.scale,corresp_tphi(slice,srow,scol,ts,trow,tcol),patchwidth/param.downscale.scale,patchheight/param.downscale.scale);
% txt=sprintf('Target: Slice %d, Row %d, Column %d',slice+1,trow,tcol);
% title(txt);

%matlabpool close
