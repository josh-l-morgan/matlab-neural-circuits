function param=photostagemap_ovupdatestagemap(param)
%This function computes an updated stage map based on the absolute
%alignment of the overview images
%By Daniel Berger for MIT-BCS Seung / Harvard Lichtman, March 2010

lparam=param.ovupdatestagemap;

for w=param.processwafers
  nrofslices=size(param.generatestagemap.linearslicepos{w},1);
  if lparam.useregionalign==1
    gabsmatrixcount=param.ovregionabsolutealign.gabsmatrixcount{w};
    gabsmatrix=param.ovregionabsolutealign.gabsmatrix{w};
  else
    gabsmatrixcount=param.ovabsolutealign.gabsmatrixcount{w};
    gabsmatrix=param.ovabsolutealign.gabsmatrix{w};
  end;
  
  in_stagemap=[param.stagemapsdir sprintf(lparam.sourcexmlfiletemplate,w)];
  out_stagemap=[param.stagemapsdir sprintf(lparam.targetxmlfiletemplate,w)];
  if lparam.overridepixelsizemu==0
    ovpixelsizenm=param.overviewimages.pixelsize_nm(w,1,1); %Use st1 sec1 to get pixel size scale factor
    ovpixelsizemu=ovpixelsizenm/1000;
  else
    ovpixelsizemu=lparam.overridepixelsizemu;
  end;
  
  txt=sprintf('Loading stage map %s ...',in_stagemap); disp(txt);
  [tree, rootname, dom]=xml_read(in_stagemap);

  numpoints=tree.NumPoints;
  txt=sprintf('Stagemap has %i points.',numpoints); disp(txt);

  refname=cell(numpoints);
%   complete=cell(numpoints);
  stagex=zeros(numpoints,1);
  stagey=zeros(numpoints,1);
%   stagez=zeros(numpoints,1);
%   stagem=zeros(numpoints,1);
%   stagetilt=zeros(numpoints,1);
  stagerot=zeros(numpoints,1);
  beamrot=zeros(numpoints,1);
%   beam_wd=zeros(numpoints,1);
%   beam_stigx=zeros(numpoints,1);
%   beam_stigy=zeros(numpoints,1);
%   detector_b=zeros(numpoints,1);
%   detector_c=zeros(numpoints,1);
%   wdcompmode=cell(numpoints);

  % Read out stage positions
  for p=1:1:numpoints %tree1.NumPoints
    tag=sprintf('Ref%i',p);
    command=sprintf('refname{p}=tree.%s.Name;',tag); eval(command);
%     command=sprintf('complete{p}=tree1.%s.Complete;',tag); eval(command);
    command=sprintf('stagex(p)=tree.%s.Stage.X;',tag); eval(command);
    command=sprintf('stagey(p)=tree.%s.Stage.Y;',tag); eval(command);
%     command=sprintf('stagez(p)=tree1.%s.Stage.Z;',tag); eval(command);
%     command=sprintf('stagem(p)=tree1.%s.Stage.M;',tag); eval(command);
%     command=sprintf('stagetilt(p)=tree1.%s.Stage.Tilt;',tag); eval(command);
    command=sprintf('stagerot(p)=tree.%s.Stage.Rot;',tag); eval(command);
    command=sprintf('beamrot(p)=tree.%s.Beam.ScanRot;',tag); eval(command);
%     command=sprintf('beam_wd(p)=tree1.%s.Beam.WD;',tag); eval(command);
%     command=sprintf('beam_stigx(p)=tree1.%s.Beam.StigX;',tag); eval(command);
%     command=sprintf('beam_stigy(p)=tree1.%s.Beam.StigY;',tag); eval(command);
%     command=sprintf('detector_b(p)=tree1.%s.Detector.B;',tag); eval(command);
%     command=sprintf('detector_c(p)=tree1.%s.Detector.C;',tag); eval(command);
%     command=sprintf('wdcompmode{p}=tree1.%s.WDcomp.ATTRIBUTE.mode;',tag); eval(command);
  end;

  %generate strip section list
%   stsec=zeros(nrofslices,2);
%   i=1;
%   nrofstrips=size(param.generatestagemap.stripslicepos{w},1);
%   for st=1:1:nrofstrips
%     nrofsections=size(param.generatestagemap.stripslicepos{w}{st},1);
%     for sec=1:1:nrofsections
%       stsec(i,1)=st; stsec(i,2)=sec;
%       i=i+1;
%     end;
%   end;
  stsec=param.alignoverviewimages.stsec{w};
  nrofremappedslices=size(stsec,1);
  
  %Update stage map
  stagexu=stagex; %in microns
  stageyu=stagey; %in microns
  stagerotu=stagerot; %in degrees
  beamrotu=beamrot; %in degrees
  
  remappedp=1;
  for p=1:1:numpoints
    if (param.generatestagemap.linearslicepos{w}(p,3)==0) %update only if it is a remapped slice (not a bad one)
      %check whether the name of the location is as expected
      expectedname=sprintf('w%02d_st%02d_sec%02d',w,stsec(remappedp,1),stsec(remappedp,2));
      actualname=refname{p};
      if(strcmp(expectedname,actualname)==0)
        txt=sprintf('WARNING: Slice name %s expected, but is %s!',expectedname,actualname); disp(txt);
      end;
      
      mtx=gabsmatrix{1,remappedp};
      beamrotu(p)=beamrotu(p)-(atan2(mtx(2,1),mtx(1,1)))*180/pi;
      stagexu(p)=stagexu(p)+mtx(1,3)*ovpixelsizemu;
      stageyu(p)=stageyu(p)-mtx(2,3)*ovpixelsizemu;
      
      remappedp=remappedp+1;
    end;
  end;
  
  
  figure(lparam.startfigure);
  plot(stagex,stagey,'*');
  hold on;
  plot(stagexu,stageyu,'r*');
  hold off;
  grid on;
  axis equal;
  title('Positions of original stagemap (blue) and updated stagemap (red)');
  
  %%%%% Construct tree for target XML file
  tree.NumPoints=numpoints;
  for pos=1:1:numpoints
    %wdcompmode=reftree.Ref1.WDcomp.ATTRIBUTE.mode; %'autofocus';
    tag=sprintf('Ref%i',pos);
    %   command=sprintf('tree.%s.Name=refname{pos};',tag); eval(command);
    %   command=sprintf('tree.%s.Complete=complete{pos};',tag); eval(command);
%     command=sprintf('tree.%s.Stage.X=stagexu(pos)+79;',tag); eval(command);
%     command=sprintf('tree.%s.Stage.Y=stageyu(pos)+114.6;',tag); eval(command);
    command=sprintf('tree.%s.Stage.X=stagexu(pos);',tag); eval(command);
    command=sprintf('tree.%s.Stage.Y=stageyu(pos);',tag); eval(command);%   command=sprintf('tree.%s.Stage.Z=stagez(pos);',tag); eval(command);
    %   command=sprintf('tree.%s.Stage.M=stagem(pos);',tag); eval(command);
    %   command=sprintf('tree.%s.Stage.Tilt=stagetilt(pos);',tag); eval(command);
    command=sprintf('tree.%s.Stage.Rot=stagerotu(pos);',tag); eval(command);
    command=sprintf('tree.%s.Beam.ScanRot=beamrotu(pos);',tag); eval(command);
    %   command=sprintf('tree.%s.Beam.WD=beam_wd(pos);',tag); eval(command);
    %   command=sprintf('tree.%s.Beam.StigX=beam_stigx(pos);',tag); eval(command);
    %   command=sprintf('tree.%s.Beam.StigY=beam_stigy(pos);',tag); eval(command);
    %   command=sprintf('tree.%s.Detector.B=detector_b(pos);',tag); eval(command);
    %   command=sprintf('tree.%s.Detector.C=detector_c(pos);',tag); eval(command);
    %command=sprintf('tree.%s.Detector.B=81.0500564575195;',tag); eval(command);
    %command=sprintf('tree.%s.Detector.C=47.9609260559082;',tag); eval(command);
    %   command=sprintf('tree.%s.WDcomp.CONTENT=[];',tag); eval(command);
    %   command=sprintf('tree.%s.WDcomp.ATTRIBUTE.mode=wdcompmode{pos};',tag); eval(command);
  end;
  %tree.MosaicSetup=tree1.MosaicSetup; %Copy all the additional parameters from the reference file
  
  txt=sprintf('Writing output XML file %s ...',out_stagemap); disp(txt);
  xml_write(out_stagemap,tree,'MosaicReferencePoints');

end;