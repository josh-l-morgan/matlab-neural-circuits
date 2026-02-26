function param=photostagemap_generatestagemap(param)
%This function is used to identify tape strips on a photographic image of a wafer.
%By Daniel Berger for MIT-BCS Seung / Harvard Lichtman, March 2010

lparam=param.generatestagemap;

txt=sprintf('Loading reference stage map %s ...',lparam.referencestagemap); disp(txt);
[reftree, refrootname, refdom]=xml_read(lparam.referencestagemap);

for w=param.processwafers
  %waferradius=param.waferradius(w);
  waferimagename=[param.normalizeddir sprintf(param.normalize.wafercroppedimagename,w)];
  manuallabelimage=[param.handlabeleddir sprintf(lparam.manuallabelimagetemplate,w)];
  targetxmlfile=[param.stagemapsdir sprintf(lparam.targetxmlfiletemplate,w)];
  nrofstrips=param.tapestrips.nrofstrips(w);
  stripxcoords=param.tapestrips.stripxcoords{w};
  stripycoords=param.tapestrips.stripycoords{w};


  if lparam.striporder==-1  %If strips are ordered right-to-left on wafer, flip strips
    stripxcoords=flipdim(stripxcoords,1);
    stripycoords=flipdim(stripycoords,1);
  end;
  
  %Display unlabeled image and overlay strip boundaries (strip 1 in blue)
  load(lparam.phototransformfile);
  oimg=imread(waferimagename);
  figure(lparam.startfigure);
  imshow(oimg);
  hold on;
  for i=1:1:nrofstrips
    x1=stripxcoords(i,1); x2=stripxcoords(i,2); x3=stripxcoords(i,3); x4=stripxcoords(i,4);
    y1=stripycoords(i,1); y2=stripycoords(i,2); y3=stripycoords(i,3); y4=stripycoords(i,4);
    if i==1
      plot([x1 x2 x4 x3 x1],[y1 y2 y4 y3 y1],'b');
    else
      plot([x1 x2 x4 x3 x1],[y1 y2 y4 y3 y1],'r');
    end;
  end;
  hold off;

  limg=imread(manuallabelimage);
  figure(lparam.startfigure+1);
  imshow(limg);
  gimg=(limg(:,:,1)==0)&(limg(:,:,2)==255)&(limg(:,:,3)==0); %All slice centers are marked in pure green
  rimg=(limg(:,:,1)==255)&(limg(:,:,2)==0)&(limg(:,:,3)==0); %All broken slice centers are marked in pure red
  
  stripnr=1:1:nrofstrips;
  cspos=cell(nrofstrips,1);
  
  
  %Create strip masks
  for i=stripnr
    x1=stripxcoords(i,1); x2=stripxcoords(i,2); x3=stripxcoords(i,3); x4=stripxcoords(i,4);
    y1=stripycoords(i,1); y2=stripycoords(i,2); y3=stripycoords(i,3); y4=stripycoords(i,4);
    xv=[x1 x2 x4 x3 x1]; yv=[y1 y2 y4 y3 y1];
    mask = poly2mask(xv,yv,size(limg,1),size(limg,2));
    gsloc=gimg.*mask;
    rsloc=rimg.*mask;
    [gL,gnum] = bwlabeln(gsloc);
    [rL,rnum] = bwlabeln(rsloc);
    txt=sprintf('Strip %i has %i good locations and %i bad locations.',i,gnum,rnum); disp(txt);
    spos=zeros(gnum+rnum,3); %y,x,color
    for p=1:1:gnum
      sbr=(gL==p);
      [yl,xl]=find(sbr);
      spos(p,1)=mean(yl);
      spos(p,2)=mean(xl);
      spos(p,3)=0; %0 is green (not broken)
    end;
    for p=1:1:rnum
      sbr=(rL==p);
      [yl,xl]=find(sbr);
      spos(gnum+p,1)=mean(yl);
      spos(gnum+p,2)=mean(xl);
      spos(gnum+p,3)=1; %1 is red (broken)
    end;
    
    spos=sortrows(spos,-1); %after this, the good and bad locations are mixed (sorted top to bottom)
    cspos{i}=spos;
    %spos
    %   figure(3+i);
    %   imagesc(stripimg);
    %   colormap gray;
  end;

  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  
  if lparam.refinemanuallabels
    % Correct label positions (of green labels only) if desired
    lrh=floor(lparam.refinetemplatesize/2);
   
    avgt=zeros(lrh*2+1,lrh*2+1);
    avgtc=zeros(lrh*2+1,lrh*2+1,3);
    avgc=0;
    
    %generate an average slice template
    for i=stripnr
     
        spos=cspos{i};
      %gnum=sum(1-spos(:,3));
      for p=1:1:size(spos,1)
         if spos(p,3)==0
          %Get current label position
          yp=spos(p,1);
          xp=spos(p,2);
          
          %Pick out region around label from unlabeled photo
         
          
          
          ir=oimg(round(yp-lrh):round(yp+lrh),round(xp-lrh):round(xp+lrh),:);        
          avgtc=avgtc+double(ir);
          figure(lparam.startfigure+2);
          subplot(1,3,1);
          imshow(ir);
          gir=mean(ir,3);
          avgt=avgt+gir; avgc=avgc+1;
          subplot(1,3,2);
          imagesc(gir); colormap gray; axis square;
          fir=bandpass2(gir,40,100);
          [A,B]=min(fir); 
          [C,D]=min(A);
          transx=D-lrh;
          transy=B(D)-lrh;
          transcorr=C;
          subplot(1,3,3);
          imagesc(fir);
          hold on;
          plot(D,B(D),'r*');
          hold off;
          %colormap jet;
          axis square;
          %input('Press Return...');
         end;
      end;
    end;
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    avgt=avgt/avgc;
    figure(lparam.startfigure+3); subplot(1,2,1);
    imagesc(avgt); colormap gray; axis square;
    title('Template patch BEFORE realign');
    
    avgtc=avgtc/avgc;
    figure(lparam.startfigure+4); imshow(uint8(avgtc));
    favgt=bandpass2(avgt,0,64);
    title('Original-color average template');
    
    ctransx=zeros(avgc,1);
    ctransy=zeros(avgc,1);
    ucspos=cspos;
    %Do template matching with average gray template avgt
    gc=1;
    for i=stripnr
      spos=cspos{i};
      uspos=ucspos{i};
      %gnum=sum(1-spos(:,3));
      for p=1:1:size(spos,1) %gnum
        if spos(p,3)==0
          %Get current label position
          yp=spos(p,1); xp=spos(p,2);
          
          %Pick out region around label from unlabeled photo
          ir=oimg(round(yp-lrh):round(yp+lrh),round(xp-lrh):round(xp+lrh),:);
          gir=mean(ir,3);
          fgir=bandpass2(gir,0,64);
          transmap=convn_fast(fgir,flipdims(favgt),'same');
          figure(lparam.startfigure+5);
          imagesc(transmap); colormap jet; axis square;
          [A,B]=max(transmap);
          [C,D]=max(A);
          transx=D-lrh-1;
          transy=B(D)-lrh-1;
          transcorr=C;
          hold on;
          plot(D,B(D),'r*');
          hold off;
          
          if (abs(transx)>lparam.refinemaxdistance)||(abs(transy)>lparam.refinemaxdistance)
            transx=0; transy=0;
          end;
          
          ctransy(gc)=transx; %D;
          ctransx(gc)=transy; %B(D);
          gc=gc+1;
          
          uspos(p,1)=uspos(p,1)+transy;
          uspos(p,2)=uspos(p,2)+transx;
          %input('Press Return...');
        end;
      end;
      ucspos{i}=uspos;
    end;
    figure(lparam.startfigure+6);
    plot(ctransx,'b');
    hold on;
    plot(ctransy,'r');
    hold off;
    grid on;
    title('Correction of marker positions based on xcorr with average');
    xlabel('Slice No.');
    ylabel('Displacement (pixels)');
      
    figure(lparam.startfigure+1);
    hold on;
    for i=stripnr
      uspos=ucspos{i};
      %gnum=sum(1-uspos(:,3));
      for p=1:1:size(uspos,1) %gnum
        plot(uspos(p,2),uspos(p,1),'r*');
      end;
    end;
    hold off;
    

    avgt2=zeros(lrh*2+1,lrh*2+1);
    avgc=0;
    for i=stripnr
      spos=ucspos{i};
      %gnum=sum(1-spos(:,3));
      for p=1:1:size(spos,1) %gnum
        if spos(p,3)==0
          %Get current label position
          yp=spos(p,1); xp=spos(p,2);
          %Pick out region around label from unlabeled photo
          ir=oimg(round(yp-lrh):round(yp+lrh),round(xp-lrh):round(xp+lrh),:);
          gir=mean(ir,3);
          avgt2=avgt2+gir; avgc=avgc+1;
        end;
      end;
    end;
    
    avgt2=avgt2/avgc;
    figure(lparam.startfigure+3); subplot(1,2,2);
    imagesc(avgt2); colormap gray; axis square;
    title('Template patch AFTER realign');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Do label position refinement again, with updated template
    favgt2=bandpass2(avgt2,0,64);
    ctransx=zeros(avgc,1);
    ctransy=zeros(avgc,1);
    ucspos=cspos;
    %Do template matching with average gray tempolate avgt
    gc=1;
    for i=stripnr
      spos=cspos{i};
      uspos=ucspos{i};
      %gnum=sum(1-spos(:,3));
      for p=1:1:size(spos,1) %gnum
        if spos(p,3)==0
          %Get current label position
          yp=spos(p,1); xp=spos(p,2);
          
          %Pick out region around label from unlabeled photo
          ir=oimg(round(yp-lrh):round(yp+lrh),round(xp-lrh):round(xp+lrh),:);
          gir=mean(ir,3);
          fgir=bandpass2(gir,0,64);
          transmap=convn_fast(fgir,flipdims(favgt2),'same');
          figure(lparam.startfigure+5);
          imagesc(transmap); colormap jet; axis square;
          [A,B]=max(transmap);
          [C,D]=max(A);
          transx=D-lrh-1;
          transy=B(D)-lrh-1;
          transcorr=C;
          hold on;
          plot(D,B(D),'r*');
          hold off;
          
          if (abs(transx)>lparam.refinemaxdistance)||(abs(transy)>lparam.refinemaxdistance)
            transx=0; transy=0;
          end;
          
          ctransy(gc)=transx; %D;
          ctransx(gc)=transy; %B(D);
          gc=gc+1;
          
          uspos(p,1)=uspos(p,1)+transy;
          uspos(p,2)=uspos(p,2)+transx;
          %input('Press Return...');
        end;
      end;
      ucspos{i}=uspos;
    end;
    figure(lparam.startfigure+7);
    plot(ctransx,'b');
    hold on;
    plot(ctransy,'r');
    hold off;
    grid on;
    title('Correction of marker positions based on xcorr with average');
    xlabel('Slice No.');
    ylabel('Displacement (pixels)');
 
    cspos=ucspos; %WRITE BACK UPDATED POSITIONS TO CSPOS FOR STAGEMAP USAGE    
    
  end; %of template-based hand-label position correction
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  nrpos=0;
  for i=stripnr
    nrpos=nrpos+size(cspos{i},1);
  end;
  txt=sprintf('This wafer has %i mapped locations.',nrpos); disp(txt);
  
  lcsname=cell(nrpos,1);
  lcspos=zeros(nrpos,3);
  
  figure(lparam.startfigure+1);
  hold on;
  lpos=1;
  for st=stripnr
    spos=cspos{st}; %this is used for the stagemap
    for pos=1:1:size(spos,1)
      if spos(pos,3)==0 %not-broken location
        name=sprintf('w%02i_st%02i_sec%02i',w,st,pos);
      else
        name=sprintf('w%02i_st%02i_sec%02i_b',w,st,pos);
      end;
      lcsname{lpos}=name; %linear list of names
      %lcspos(lpos,1)=spos(pos,1); lcspos(lpos,2)=spos(pos,2); %linear list of positions and slice state (0=good, 1=bad)
      lcspos(lpos,:)=spos(pos,:);
      pname=['  ',strrep(name,'_','\_')];
      text(spos(pos,2),spos(pos,1),pname);
      lpos=lpos+1;
    end;
  end;
  hold off;
  
  param.generatestagemap.stripslicepos{w}=cspos;
  param.generatestagemap.linearslicepos{w}=lcspos;
  param.generatestagemap.linearslicename{w}=lcsname;
  
  
  % cstripimg=stripimg([floor(min(yv)):floor(max(yv)+1)],[floor(min(xv)):floor(max(xv)+1)]);
  % figure(17);
  % imagesc(cstripimg);
  % colormap gray;
  
  
  %%%%% Construct tree for XML file
  tree.NumPoints=nrpos;
  for pos=1:1:nrpos
    name=lcsname{pos};
    iscomplete=reftree.Ref1.Complete; %'false';
    stagex=lcspos(pos,2);
    stagey=lcspos(pos,1);
    
    %Apply image->stage transformation
    dpos=T*[stagey stagex 1]';
    stagey=dpos(1); stagex=dpos(2);
    
    stagez= 7299.9; %reftree.Ref1.Stage.Z; %727.138;
    stagem=reftree.Ref1.Stage.M; %0;
    stagetilt=0; %reftree.Ref1.Stage.Tilt; %-1;
    stagerot=reftree.Ref1.Stage.Rot; %0; %Here the stage rotation of the image->stage mapping is used
    beamscanrot=reftree.Ref1.Beam.ScanRot; %0;
    beamwd=reftree.Ref1.Beam.WD; %0.0075;
    beamstigx=reftree.Ref1.Beam.StigX; %-0.32;
    beamstigy=reftree.Ref1.Beam.StigY; %-0.75;
    detectorb=reftree.Ref1.Detector.B; %69.2;
    detectorc=reftree.Ref1.Detector.C; %38.05;
    wdcompmode=reftree.Ref1.WDcomp.ATTRIBUTE.mode; %'autofocus';
    
    tag=sprintf('Ref%i',pos);
    command=sprintf('tree.%s.Name=name;',tag); eval(command);
    command=sprintf('tree.%s.Complete=iscomplete;',tag); eval(command);
    
    command=sprintf('tree.%s.Stage.X=stagex;',tag); eval(command);
    command=sprintf('tree.%s.Stage.Y=stagey;',tag); eval(command);
    command=sprintf('tree.%s.Stage.Z=stagez;',tag); eval(command);
    command=sprintf('tree.%s.Stage.M=stagem;',tag); eval(command);
    command=sprintf('tree.%s.Stage.Tilt=stagetilt;',tag); eval(command);
    command=sprintf('tree.%s.Stage.Rot=stagerot;',tag); eval(command);
    command=sprintf('tree.%s.Beam.ScanRot=beamscanrot;',tag); eval(command);
    command=sprintf('tree.%s.Beam.WD=beamwd;',tag); eval(command);
    command=sprintf('tree.%s.Beam.StigX=beamstigx;',tag); eval(command);
    command=sprintf('tree.%s.Beam.StigY=beamstigy;',tag); eval(command);
    command=sprintf('tree.%s.Detector.B=detectorb;',tag); eval(command);
    command=sprintf('tree.%s.Detector.C=detectorc;',tag); eval(command);
    command=sprintf('tree.%s.WDcomp.CONTENT=[];',tag); eval(command);
    command=sprintf('tree.%s.WDcomp.ATTRIBUTE.mode=wdcompmode;',tag); eval(command);
  end;
  tree.MosaicSetup=reftree.MosaicSetup; %Copy all the additional parameters from the reference file
  
  txt=sprintf('Writing XML file %s ...',targetxmlfile); disp(txt);
  xml_write(targetxmlfile,tree,'MosaicReferencePoints');
  
  if (lparam.stopeach)&&(w<param.processwafers(end))
    input('Press Return...');
  end;
end;

%Compute and store overall slice count
param.allwafergoodslicecount=0;
param.allwaferbadslicecount=0;
for w=1:1:size(param.generatestagemap.linearslicepos,2)
  a=param.generatestagemap.linearslicepos{w};
  for s=1:1:size(a,1)
    if a(s,3)==0
      param.allwafergoodslicecount=param.allwafergoodslicecount+1;
    else
      param.allwaferbadslicecount=param.allwaferbadslicecount+1;
    end;
  end;
end;
txt=sprintf('The %d mapped wafers contain %d good and %d bad slices overall.',size(param.generatestagemap.linearslicepos,2),param.allwafergoodslicecount,param.allwaferbadslicecount);
disp(txt);