function param=pipeline_fititerativeaffine(param)
%This script computes a global affine fit of a volume of tiles, based
%on previously computed local correspondences.
%For use with the pipeline.
%By Daniel Berger for MIT-BCS Seung, June 22 2009

lparam=param.fititerativeaffine;
imagewidth=param.rawsize(2);
imageheight=param.rawsize(1);
fakt=1;
fimagewidth=imagewidth*fakt;
fimageheight=imageheight*fakt;
gridwidth=param.predictcorr.gridwidth;
gridheight=param.predictcorr.gridheight;

corresp_corr_sx=param.computelocalcorr.corresp_csx;
corresp_corr_sy=param.computelocalcorr.corresp_csy;
corresp_corr_tx=param.computelocalcorr.corresp_ctx;
corresp_corr_ty=param.computelocalcorr.corresp_cty;
corresp_corr_max=param.computelocalcorr.corresp_cmax;
corresp_corr_vargt075=param.computelocalcorr.corresp_cvar075;

disp('Initializing affine alignment with rigid alignment...');
affine=param.rigid.absrigidmatrix;
%The rigid matrix is in downscaled coordinates; scale up translations
for slice=1:1:param.nrofslices
  for row=1:1:param.nrofrows
    for column=1:1:param.nrofcolumns
      affine(slice,row,column,1,3)=affine(slice,row,column,1,3)*param.downscale.scale;
      affine(slice,row,column,2,3)=affine(slice,row,column,2,3)*param.downscale.scale;
    end;
  end;
end;

disp('Computing all global coordinates of corresponding points...');
%Compute global coordinates of all corresponding pairs, using the
%pre-initialized affine transformations
global cp_gsx cp_gsy cp_gtx cp_gty;
cp_gsx=zeros(param.nrofslices,param.nrofrows,param.nrofcolumns,3,param.nrofrows,param.nrofcolumns,gridheight,gridwidth);
cp_gsy=zeros(param.nrofslices,param.nrofrows,param.nrofcolumns,3,param.nrofrows,param.nrofcolumns,gridheight,gridwidth);
cp_gtx=zeros(param.nrofslices,param.nrofrows,param.nrofcolumns,3,param.nrofrows,param.nrofcolumns,gridheight,gridwidth);
cp_gty=zeros(param.nrofslices,param.nrofrows,param.nrofcolumns,3,param.nrofrows,param.nrofcolumns,gridheight,gridwidth);

for s=1:1:param.nrofslices %slice
  for sy=1:1:param.nrofrows %source row
    for sx=1:1:param.nrofcolumns %source column
      sA=squeeze(affine(s,sy,sx,:,:)); %affine transform for source tile is constant for all further loops
      for ts=1:1:3 %target slice (1=prev, 2=same, 3=next)
        if (s+ts-2>0) && (s+ts-2<=param.nrofslices) %only if target slice exists
          for ty=1:1:param.nrofrows %target row
            for tx=1:1:param.nrofcolumns %target column
              tA=squeeze(affine(s+ts-2,ty,tx,:,:)); %affine transformation for target tile
              for gy=1:1:gridheight %grid row
                for gx=1:1:gridwidth %grid column
                  if (corresp_corr_vargt075(s,sy,sx,ts,ty,tx,gy,gx)<lparam.vargt075threshold) && (corresp_corr_vargt075(s,sy,sx,ts,ty,tx,gy,gx)>0)  %only use points which are robust
                    xc=corresp_corr_sx(s,sy,sx,ts,ty,tx,gy,gx)*fakt;
                    yc=corresp_corr_sy(s,sy,sx,ts,ty,tx,gy,gx)*fakt;
                    if (xc~=0) && (yc~=0) %olny actually existing corresponding points are transformed
%                       gc=sA*[xc, yc, 1]';
%                       cp_gsx(s,sy,sx,ts,ty,tx,gy,gx)=gc(1);
%                       cp_gsy(s,sy,sx,ts,ty,tx,gy,gx)=gc(2);
                      %CAUTION - We have to turn around each tile image center.
                      gc=sA*[xc-fimagewidth/2, yc-fimageheight/2, 1]';
                      cp_gsx(s,sy,sx,ts,ty,tx,gy,gx)=gc(1)+fimagewidth/2;
                      cp_gsy(s,sy,sx,ts,ty,tx,gy,gx)=gc(2)+fimageheight/2;
                    end;
                    xc=corresp_corr_tx(s,sy,sx,ts,ty,tx,gy,gx)*fakt;
                    yc=corresp_corr_ty(s,sy,sx,ts,ty,tx,gy,gx)*fakt;
                    if (xc~=0) && (yc~=0)
                      gc=tA*[xc-fimagewidth/2, yc-fimageheight/2, 1]';
                      cp_gtx(s,sy,sx,ts,ty,tx,gy,gx)=gc(1)+fimagewidth/2;
                      cp_gty(s,sy,sx,ts,ty,tx,gy,gx)=gc(2)+fimageheight/2;
                    end;
                  end;
                end;
              end;
            end;
          end;
        end;
      end;
    end;
  end;
end;

%showprealignment=0;

if lparam.showprealignment==1
  %plot slices and global in-slice corresponding point pairs
  figure(1);
  for s=1:1:param.nrofslices
    for row=1:1:param.nrofrows
      for column=1:1:param.nrofcolumns
        A=squeeze(affine(s,row,column,:,:));
        %coord=A*[[1 fimagewidth fimagewidth 1]; [1 1 fimageheight fimageheight]; [1 1 1 1]];
        coord=A*[[-fimagewidth/2 fimagewidth/2 fimagewidth/2 -fimagewidth/2]; [-fimageheight/2 -fimageheight/2 fimageheight/2 fimageheight/2]; [1 1 1 1]];
        drawpoly(coord(1,:)'+fimagewidth/2,coord(2,:)'+fimageheight/2,'k');
        hold on;
        if s>1
          A=squeeze(affine(s-1,row,column,:,:));
          %coord=A*[[1 fimagewidth fimagewidth 1]; [1 1 fimageheight fimageheight]; [1 1 1 1]];
          coord=A*[[-fimagewidth/2 fimagewidth/2 fimagewidth/2 -fimagewidth/2]; [-fimageheight/2 -fimageheight/2 fimageheight/2 fimageheight/2]; [1 1 1 1]];
          drawpoly(coord(1,:)'+fimagewidth/2,coord(2,:)'+fimageheight/2,'b--');
          hold on;
        end;
        if s<param.nrofslices
          A=squeeze(affine(s+1,row,column,:,:));
          %coord=A*[[1 fimagewidth fimagewidth 1]; [1 1 fimageheight fimageheight]; [1 1 1 1]];
          coord=A*[[-fimagewidth/2 fimagewidth/2 fimagewidth/2 -fimagewidth/2]; [-fimageheight/2 -fimageheight/2 fimageheight/2 fimageheight/2]; [1 1 1 1]];
          drawpoly(coord(1,:)'+fimagewidth/2,coord(2,:)'+fimageheight/2,'g--');
          hold on;
        end;
        
        %same-slice
        for ty=1:1:param.nrofrows
          for tx=1:1:param.nrofcolumns
            if (ty~=row) || (tx~=column)
              for gy=1:1:gridheight
                for gx=1:1:gridwidth
                  if cp_gtx(s,row,column,2,ty,tx,gy,gx)~=0
                    gsx=cp_gsx(s,row,column,2,ty,tx,gy,gx);
                    gsy=cp_gsy(s,row,column,2,ty,tx,gy,gx);
                    gtx=cp_gtx(s,row,column,2,ty,tx,gy,gx);
                    gty=cp_gty(s,row,column,2,ty,tx,gy,gx);
                    if corresp_corr_vargt075(s,row,column,2,ty,tx,gy,gx)<lparam.vargt075threshold
                      plot(gsx,gsy,'k.');
                      plot(gtx,gty,'k*');
                      plot([gsx gtx],[gsy gty],'k');
                    else %this will never be shown because global coordinates are only computed for points below vargt075 threshold
                      plot(gsx,gsy,'r.');
                      plot(gtx,gty,'r*');
                      plot([gsx gtx],[gsy gty],'r');
                    end;
                  end;
                end;
              end;
            end;
          end;
        end;
        
        %previous slice
        for ty=1:1:param.nrofrows
          for tx=1:1:param.nrofcolumns
            %if (ty~=row) || (tx~=column)
            for gy=1:1:gridheight
              for gx=1:1:gridwidth
                if cp_gtx(s,row,column,1,ty,tx,gy,gx)~=0
                  gsx=cp_gsx(s,row,column,1,ty,tx,gy,gx);
                  gsy=cp_gsy(s,row,column,1,ty,tx,gy,gx);
                  gtx=cp_gtx(s,row,column,1,ty,tx,gy,gx);
                  gty=cp_gty(s,row,column,1,ty,tx,gy,gx);
                  if corresp_corr_vargt075(s,row,column,1,ty,tx,gy,gx)<lparam.vargt075threshold
                    plot(gsx,gsy,'b.');
                    plot(gtx,gty,'bx');
                    plot([gsx gtx],[gsy gty],'b');
                  else
                    plot(gsx,gsy,'r.');
                    plot(gtx,gty,'rx');
                    plot([gsx gtx],[gsy gty],'r');
                  end;
                end;
              end;
              %end;
            end;
          end;
        end;
        
        %next slice
        for ty=1:1:param.nrofrows
          for tx=1:1:param.nrofcolumns
            %if (ty~=row) || (tx~=column)
            for gy=1:1:gridheight
              for gx=1:1:gridwidth
                if cp_gtx(s,row,column,3,ty,tx,gy,gx)~=0
                  gsx=cp_gsx(s,row,column,3,ty,tx,gy,gx);
                  gsy=cp_gsy(s,row,column,3,ty,tx,gy,gx);
                  gtx=cp_gtx(s,row,column,3,ty,tx,gy,gx);
                  gty=cp_gty(s,row,column,3,ty,tx,gy,gx);
                  if corresp_corr_vargt075(s,row,column,3,ty,tx,gy,gx)<lparam.vargt075threshold
                    plot(gsx,gsy,'g.');
                    plot(gtx,gty,'go');
                    plot([gsx gtx],[gsy gty],'g');
                  else
                    plot(gsx,gsy,'r.');
                    plot(gtx,gty,'ro');
                    plot([gsx gtx],[gsy gty],'r');
                  end;
                end;
              end;
              %end;
            end;
          end;
        end;
      end;
    end;
    hold off;
    txt=sprintf('Slice No. %d ',s);
    title(txt);
    axis equal; %square;
    %axis([-4000*fakt 14000*fakt -4000*fakt 14000*fakt]);
    grid on;
    %if s==305
      input('press RETURN');
    %end;
  end;
  
end;


disp('Refine affine alignment iteratively...');
%nrofiterations=25000; %this is now passed as lparam.nrofiterations

affparams=zeros(param.nrofslices*param.nrofrows*param.nrofcolumns,6);
meanaffparams=zeros(lparam.nrofiterations,6);
meandist=zeros(lparam.nrofiterations,1);

%perm stores the current random permutation
if lparam.fixfirst==1
  perm=zeros(param.nrofslices*param.nrofrows*param.nrofcolumns-1,4); %-1 if you want to optimize all except tile(1,1,1)
else
  perm=zeros(param.nrofslices*param.nrofrows*param.nrofcolumns,4); %-1 if you want to optimize all except tile(1,1,1)
end;
  %Iteratively refine affine transformations in random order
for iteration=1:1:lparam.nrofiterations
  txt=sprintf('Iteration %d/%d...',iteration,lparam.nrofiterations);
  disp(txt);

  %create a randomly ordered list of tiles
  p=1;
  for slice=1:1:param.nrofslices
    for row=1:1:param.nrofrows
      for column=1:1:param.nrofcolumns
        if (lparam.fixfirst==0) || (slice~=1) || (row~=1) || (column~=1) %optimize all except tile(1,1,1)
          perm(p,1)=slice;
          perm(p,2)=row;
          perm(p,3)=column;
          perm(p,4)=random('unif',0,1000);
          p=p+1;
        end;
      end;
    end;
  end;
  perm=sortrows(perm,4);

    
  cmeandist=0;
  cmeandistc=0;
  
  %go through list and optimize each affine transformation with respect to
  %its immediate neighbors (which are kept fixed)
  for i=1:1:size(perm,1)
    %   if floor(i/10)==i/10
    %     disp(i);
    %   end;
    as=perm(i,1);
    arow=perm(i,2);
    acolumn=perm(i,3);
    
%      if (iteration==5) && (as==2) && (arow==2) && (acolumn==2)
%        input('press return');
%      end;
    
    %needs global cp_gsx cp_gsy cp_gtx cp_gty;
    [cpsx,cpsy,cptx,cpty]=computeaffine_iterative_getcplist_nooutliers(param.nrofslices, param.nrofrows, param.nrofcolumns, gridwidth,gridheight, as,arow,acolumn);
    
    cmeandist=cmeandist+sum(((cpty-cpsy).*(cpty-cpsy)+(cptx-cpsx).*(cptx-cpsx)));
    cmeandistc=cmeandistc+size(cpsx,2);
    
%     A=getAffine(cpsx,cpsy,cptx,cpty); %A is the GLOBAL transform (inverse of local)
%     %A=(A+[[1 0 0]; [0 1 0]; [0 0 1]])/2; %only go HALF the way (this prevents oscillation; is that the correct solution..?)
%     sA=squeeze(affine(as,arow,acolumn,:,:));

    A=getAffine(cpsx-fimagewidth/2,cpsy-fimageheight/2,cptx-fimagewidth/2,cpty-fimageheight/2); %A is the GLOBAL transform (inverse of local)
    %A(1,3)=A(1,3)+fimagewidth/2;
    %A(2,3)=A(2,3)+fimageheight/2;
    %A=(A+[[1 0 0]; [0 1 0]; [0 0 1]])/2; %only go HALF the way (this prevents oscillation; is that the correct solution..?)
    sA=squeeze(affine(as,arow,acolumn,:,:));

    tA=A*sA; %sA*A;
    [affparams(i,1),affparams(i,2),affparams(i,3),affparams(i,4),affparams(i,5),affparams(i,6)]=getfromaffinematrix(A);
    
    if lparam.tweakaffine==1
      [scalex,scaley,shearx,rotang,transx,transy]=getfromaffinematrix(tA);
      shearx=shearx*(1-lparam.tweakshear); %0.95; 
      %isotropy: scalex and scaley should be similar
      avgscale=(scalex+scaley)/2;
      aspectratio=scalex/scaley;
      regavgscale=avgscale*(1-lparam.tweakscale) + 1*lparam.tweakscale;
      regaspectratio=aspectratio*(1-lparam.tweakiso) + 1*lparam.tweakiso;
      scaley=(2*regavgscale)/(regaspectratio+1);
      scalex=regaspectratio*scaley;
      
%       scalexiso=scalex*(1-lparam.tweakiso) + avgscale*lparam.tweakiso;
%       scaleyiso=scaley*(1-lparam.tweakiso) + avgscale*lparam.tweakiso;
%       scalex=scalexiso*(1-lparam.tweakscale) + 1*lparam.tweakscale; %0.95+0.05; 
%       scaley=scaleyiso*(1-lparam.tweakscale) + 1*lparam.tweakscale; %*0.95+0.05;
%       avgscale=(avgscale*(1-lparam.tweakscale))+(1*lparam.tweakscale); %relax average scale towards 1 by tweakscale
%       scalex=(scalex*(1-lparam.tweakiso))+(avgscale*lparam.tweakiso); %relax scalex towards average scale by tweakiso
%       scaley=(scaley*(1-lparam.tweakiso))+(avgscale*lparam.tweakiso); %relax scaley towards average scale by tweakiso
      %scalex=scalex*lparam.tweakscaleval+(1-lparam.tweakscaleval); %0.95+0.05; 
      %scaley=scaley*lparam.tweakscaleval+(1-lparam.tweakscaleval); %*0.95+0.05;
      tA=makeaffinematrix(scalex,scaley,shearx,0,rotang,transx,transy);
    end;
    
%     
%     %if (as==20) && (arow==1) && (acolumn==2)
%       %Visualize change of affine mapping
%       %figure(floor(10+iteration/100));
%       figure(10+as*4+arow*2+acolumn);
%       coord=sA*[[1 imagewidth imagewidth 1]; [1 1 imageheight imageheight]; [1 1 1 1]];
%       drawpoly(coord(1,:)',coord(2,:)','b');
%       hold on;
%       coord=tA*[[1 imagewidth imagewidth 1]; [1 1 imageheight imageheight]; [1 1 1 1]];
%       drawpoly(coord(1,:)',coord(2,:)','k');
%       plot([cpsx; cptx],[cpsy; cpty],'b');
%       plot(cptx,cpty,'b.');
%       %Ap=A*[cptx; cpty; ones(1,size(cptx,2))];
%       %plot(cpsx,cpsy,'ko');
%       %plot([cpsx; Ap(1,:)],[cpsy; Ap(2,:)],'k');
%       %plot(Ap(1,:),Ap(2,:),'k.');
%     %end;
    
    %update affine transformation
    affine(as,arow,acolumn,:,:)=tA;
    
    %update corresponding-point-locations
 
    
    %previous-slice tiles
    s=as-1;
    if (s>0) && (s<=param.nrofslices)
      for sy=arow-1:arow+1 %source row
        if (sy>0) && (sy<=param.nrofrows)
          for sx=acolumn-1:acolumn+1 %source column
            if (sx>0) && (sx<=param.nrofcolumns)
              for gy=1:1:gridheight %grid row
                for gx=1:1:gridwidth %grid column
                  %previous-slice tiles -> updated tile
                  if (corresp_corr_vargt075(s,sy,sx,3,arow,acolumn,gy,gx)<lparam.vargt075threshold) && (corresp_corr_vargt075(s,sy,sx,3,arow,acolumn,gy,gx)>0)  %only use points which are robust
                    xc=corresp_corr_tx(s,sy,sx,3,arow,acolumn,gy,gx)*fakt;
                    yc=corresp_corr_ty(s,sy,sx,3,arow,acolumn,gy,gx)*fakt;
                    if (xc~=0) && (yc~=0)
                      %ps=(cp_gtx(s,sy,sx,3,arow,acolumn,gy,gx)-cp_gsx(s,sy,sx,3,arow,acolumn,gy,gx))
                      %predistsum=predistsum+ps;
                      gc=tA*[xc-fimagewidth/2, yc-fimageheight/2, 1]';
                      cp_gtx(s,sy,sx,3,arow,acolumn,gy,gx)=gc(1)+fimagewidth/2;
                      cp_gty(s,sy,sx,3,arow,acolumn,gy,gx)=gc(2)+fimageheight/2;
                    end;
                  end;
                  %updated tile -> previous-slice tiles
                  if (corresp_corr_vargt075(as,arow,acolumn,1,sy,sx,gy,gx)<lparam.vargt075threshold) && (corresp_corr_vargt075(as,arow,acolumn,1,sy,sx,gy,gx)>0)  %only use points which are robust
                    xc=corresp_corr_sx(as,arow,acolumn,1,sy,sx,gy,gx)*fakt;
                    yc=corresp_corr_sy(as,arow,acolumn,1,sy,sx,gy,gx)*fakt;
                    if (xc~=0) && (yc~=0)
                      gc=tA*[xc-fimagewidth/2, yc-fimageheight/2, 1]';
                      cp_gsx(as,arow,acolumn,1,sy,sx,gy,gx)=gc(1)+fimagewidth/2;
                      cp_gsy(as,arow,acolumn,1,sy,sx,gy,gx)=gc(2)+fimageheight/2;
                    end;
                  end;
                end;
              end;
            end;
          end;
        end;
      end;
    end;
    
%     if (iteration==32) && (as==1) && (arow==2) && (acolumn==2)
%       input('press key 2');
%     end;
    
    %same-slice tiles
    s=as;
    for sy=arow-1:arow+1 %source row
      if (sy>0) && (sy<=param.nrofrows)
        for sx=acolumn-1:acolumn+1 %source column
          if (sx>0) && (sx<=param.nrofcolumns) && ((sy~=arow) || (sx~=acolumn)) %all except the same tile
            for gy=1:1:gridheight %grid row
              for gx=1:1:gridwidth %grid column
                %same-slice tiles -> updated tile
%                 if (iteration==32) && (as==1) && (arow==2) && (acolumn==2) && (sy==1) && (sx==2) && (gy==1) && (gx==8)
%                   input('press key 2');
%                 end;
                
                if (corresp_corr_vargt075(s,sy,sx,2,arow,acolumn,gy,gx)<lparam.vargt075threshold) && (corresp_corr_vargt075(s,sy,sx,2,arow,acolumn,gy,gx)>0)  %only use points which are robust
                  xc=corresp_corr_tx(s,sy,sx,2,arow,acolumn,gy,gx)*fakt;
                  yc=corresp_corr_ty(s,sy,sx,2,arow,acolumn,gy,gx)*fakt;
                  if (xc~=0) && (yc~=0)
                    gc=tA*[xc-fimagewidth/2, yc-fimageheight/2, 1]';
                    cp_gtx(s,sy,sx,2,arow,acolumn,gy,gx)=gc(1)+fimagewidth/2;
                    cp_gty(s,sy,sx,2,arow,acolumn,gy,gx)=gc(2)+fimageheight/2;
                  end;
                end;
                %updated tile -> same-slice tiles
                if (corresp_corr_vargt075(as,arow,acolumn,2,sy,sx,gy,gx)<lparam.vargt075threshold) && (corresp_corr_vargt075(as,arow,acolumn,2,sy,sx,gy,gx)>0)  %only use points which are robust
                  xc=corresp_corr_sx(as,arow,acolumn,2,sy,sx,gy,gx)*fakt;
                  yc=corresp_corr_sy(as,arow,acolumn,2,sy,sx,gy,gx)*fakt;
                  if (xc~=0) && (yc~=0)
                    gc=tA*[xc-fimagewidth/2, yc-fimageheight/2, 1]';
                    cp_gsx(as,arow,acolumn,2,sy,sx,gy,gx)=gc(1)+fimagewidth/2;
                    cp_gsy(as,arow,acolumn,2,sy,sx,gy,gx)=gc(2)+fimageheight/2;
                  end;
                end;
              end;
            end;
          end;
        end;
      end;
    end;
    
    %next-slice tiles
    s=as+1;
    if (s>0) && (s<=param.nrofslices)
      for sy=arow-1:arow+1 %source row
        if (sy>0) && (sy<=param.nrofrows)
          for sx=acolumn-1:acolumn+1 %source column
            if (sx>0) && (sx<=param.nrofcolumns)
              for gy=1:1:gridheight %grid row
                for gx=1:1:gridwidth %grid column
                  %next-slice tiles -> updated tile
                  if (corresp_corr_vargt075(s,sy,sx,1,arow,acolumn,gy,gx)<lparam.vargt075threshold) && (corresp_corr_vargt075(s,sy,sx,1,arow,acolumn,gy,gx)>0)  %only use points which are robust
                    xc=corresp_corr_tx(s,sy,sx,1,arow,acolumn,gy,gx)*fakt;
                    yc=corresp_corr_ty(s,sy,sx,1,arow,acolumn,gy,gx)*fakt;
                    if (xc~=0) && (yc~=0)
                      gc=tA*[xc-fimagewidth/2, yc-fimageheight/2, 1]';
                      cp_gtx(s,sy,sx,1,arow,acolumn,gy,gx)=gc(1)+fimagewidth/2;
                      cp_gty(s,sy,sx,1,arow,acolumn,gy,gx)=gc(2)+fimageheight/2;
                    end;
                  end;
                  %updated tile -> next-slice tiles
                  if (corresp_corr_vargt075(as,arow,acolumn,3,sy,sx,gy,gx)<lparam.vargt075threshold) && (corresp_corr_vargt075(as,arow,acolumn,3,sy,sx,gy,gx)>0)  %only use points which are robust
                    xc=corresp_corr_sx(as,arow,acolumn,3,sy,sx,gy,gx)*fakt;
                    yc=corresp_corr_sy(as,arow,acolumn,3,sy,sx,gy,gx)*fakt;
                    if (xc~=0) && (yc~=0)
                      gc=tA*[xc-fimagewidth/2, yc-fimageheight/2, 1]';
                      cp_gsx(as,arow,acolumn,3,sy,sx,gy,gx)=gc(1)+fimagewidth/2;
                      cp_gsy(as,arow,acolumn,3,sy,sx,gy,gx)=gc(2)+fimageheight/2;
                    end;
                  end;
                end;
              end;
            end;
          end;
        end;
      end;
    end;
    
    [ncpsx,ncpsy,ncptx,ncpty]=computeaffine_iterative_getcplist_nooutliers(param.nrofslices,param.nrofrows,param.nrofcolumns, gridwidth,gridheight, as,arow,acolumn);
        
    predistsum=sum((cptx-cpsx).*(cptx-cpsx) + (cpty-cpsy).*(cpty-cpsy));
    postdistsum=sum((ncptx-ncpsx).*(ncptx-ncpsx) + (ncpty-ncpsy).*(ncpty-ncpsy));
    %postdistsum-predistsum

    
    
    %   figure(3);
    %   s=as;
    %   for row=1:1:nrofrows
    %     for column=1:1:nrofcolumns
    %       A=squeeze(affine(s,row,column,:,:));
    %       coord=A*[[1 imagewidth imagewidth 1]; [1 1 imageheight imageheight]; [1 1 1 1]];
    %       drawpoly(coord(1,:)',coord(2,:)','k');
    %       hold on;
    % %       if s>1
    % %         A=squeeze(affine(s-1,row,column,:,:));
    % %         coord=A*[[1 imagewidth imagewidth 1]; [1 1 imageheight imageheight]; [1 1 1 1]];
    % %         drawpoly(coord(1,:)',coord(2,:)','b--');
    % %         hold on;
    % %       end;
    %
    %       %same-slice
    %       for ty=1:1:nrofrows
    %         for tx=1:1:nrofcolumns
    %           if (ty~=row) || (tx~=column)
    %             for gy=1:1:8
    %               for gx=1:1:8
    %                 if cp_gtx(s,row,column,2,ty,tx,gy,gx)~=0
    %                   gsx=cp_gsx(s,row,column,2,ty,tx,gy,gx);
    %                   gsy=cp_gsy(s,row,column,2,ty,tx,gy,gx);
    %                   gtx=cp_gtx(s,row,column,2,ty,tx,gy,gx);
    %                   gty=cp_gty(s,row,column,2,ty,tx,gy,gx);
    %                   if lin_corresp_corr_vargt075(s,row,column,2,ty,tx,gy,gx)<threshold
    %                     plot(gsx,gsy,'k.');
    %                     plot(gtx,gty,'k*');
    %                     plot([gsx gtx],[gsy gty],'k');
    %                   else
    %                     plot(gsx,gsy,'r.');
    %                     plot(gtx,gty,'r*');
    %                     plot([gsx gtx],[gsy gty],'r');
    %                   end;
    %                 end;
    %               end;
    %             end;
    %           end;
    %         end;
    %       end;
    %
    %       %previous slice
    %       for ty=1:1:nrofrows
    %         for tx=1:1:nrofcolumns
    %           %if (ty~=row) || (tx~=column)
    %           for gy=1:1:8
    %             for gx=1:1:8
    %               if cp_gtx(s,row,column,1,ty,tx,gy,gx)~=0
    %                 gsx=cp_gsx(s,row,column,1,ty,tx,gy,gx);
    %                 gsy=cp_gsy(s,row,column,1,ty,tx,gy,gx);
    %                 gtx=cp_gtx(s,row,column,1,ty,tx,gy,gx);
    %                 gty=cp_gty(s,row,column,1,ty,tx,gy,gx);
    %                 if lin_corresp_corr_vargt075(s,row,column,1,ty,tx,gy,gx)<threshold
    %                   plot(gsx,gsy,'b.');
    %                   plot(gtx,gty,'b*');
    %                   plot([gsx gtx],[gsy gty],'b');
    %                 else
    %                   plot(gsx,gsy,'r.');
    %                   plot(gtx,gty,'r*');
    %                   plot([gsx gtx],[gsy gty],'r');
    %                 end;
    %               end;
    %             end;
    %             %end;
    %           end;
    %         end;
    %       end;
    %
    %       %next slice
    %       for ty=1:1:nrofrows
    %         for tx=1:1:nrofcolumns
    %           %if (ty~=row) || (tx~=column)
    %           for gy=1:1:8
    %             for gx=1:1:8
    %               if cp_gtx(s,row,column,3,ty,tx,gy,gx)~=0
    %                 gsx=cp_gsx(s,row,column,3,ty,tx,gy,gx);
    %                 gsy=cp_gsy(s,row,column,3,ty,tx,gy,gx);
    %                 gtx=cp_gtx(s,row,column,3,ty,tx,gy,gx);
    %                 gty=cp_gty(s,row,column,3,ty,tx,gy,gx);
    %                 if lin_corresp_corr_vargt075(s,row,column,3,ty,tx,gy,gx)<threshold
    %                   plot(gsx,gsy,'g.');
    %                   plot(gtx,gty,'go');
    %                   plot([gsx gtx],[gsy gty],'g');
    %                 else
    %                   plot(gsx,gsy,'r.');
    %                   plot(gtx,gty,'ro');
    %                   plot([gsx gtx],[gsy gty],'r');
    %                 end;
    %               end;
    %             end;
    %             %end;
    %           end;
    %         end;
    %       end;
    %
    %     end;
    %   end;
    %   hold off;
    %   txt=sprintf('Slice %d, Row %d, Column %d BEFORE',as,arow,acolumn);
    %   title(txt);
    %input('Press RETURN');
  end;
  
  %Visualize progress
  meanaffparams(iteration,:)=mean(affparams,1);
  figure(floor(iteration/100)+30);
  hist(affparams,50);
  txt=sprintf('Histogram of affine transform parameters for iteration %d',iteration);
  title(txt);
  
  figure(5);
  plot(meanaffparams(1:iteration,:));
  grid on;
  txt=sprintf('Mean affine transform parameters for iterations 1 - %d',iteration);
  title(txt);
  
  meandist(iteration)=cmeandist/cmeandistc;
  figure(6);
  semilogy(meandist(1:iteration));
  grid on;
  txt=sprintf('Mean corresponding points distance for iterations 1 - %d',iteration);
  title(txt);
  
%   if (iteration/100)==floor(iteration/100) %save every 100th iteration
%     save('affine_dunoit_fixed.mat','affine','meanaffparams','slicenr','wafernr');
%   end;  
  param.fititerativeaffine.affine=affine;
end;


if lparam.postprocessing==1 %if this is 1, postprocessing is performed
  %read out average scale and rotation angle
  avgscale=0; avgrot=0; count=0; 
  for slice=1:1:param.nrofslices
    for row=1:1:param.nrofrows
      for column=1:1:param.nrofcolumns
        A=squeeze(param.fititerativeaffine.affine(slice,row,column,:,:));
        [scalex,scaley,shearx,rotang,transx,transy]=getfromaffinematrix(A);
        avgscale=avgscale+(scalex+scaley)/2;
        avgrot=avgrot+rotang;
        count=count+1;
      end;
    end;
  end;
  avgscale=avgscale/count;
  avgrot=avgrot/count;
  
  if (lparam.post_scale==0)
    avgscale=1;
  end;
  
  if (lparam.post_rotate==0)
    avgrot=0;
  end;
  
  rigidmtx=[[cos(-avgrot)/avgscale -sin(-avgrot)/avgscale 0]; [sin(-avgrot)/avgscale cos(-avgrot)/avgscale 0]; [0 0 1]];
  
  %correct affine matrices
  for slice=1:1:param.nrofslices
    for row=1:1:param.nrofrows
      for column=1:1:param.nrofcolumns
        A=squeeze(param.fititerativeaffine.affine(slice,row,column,:,:));
%         [scalex,scaley,shearx,rotang,transx,transy]=getfromaffinematrix(A);
%         if lparam.post_scale==1     %normalize overall scaling
%           scalex=scalex/avgscale;
%           scaley=scaley/avgscale;
%         end;
%         if lparam.post_rotate==1    %normalize overall rotation
%           rotang=rotang-avgrot;
%         end;
%         param.fititerativeaffine.affine(slice,row,column,:,:)=makeaffinematrix(scalex,scaley,shearx,0,rotang,transx,transy);
        param.fititerativeaffine.affine(slice,row,column,:,:)=rigidmtx*A;
      end;
    end;
  end;
end;

