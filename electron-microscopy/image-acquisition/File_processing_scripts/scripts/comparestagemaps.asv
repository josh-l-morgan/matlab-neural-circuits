function diffpos=comparestagemaps(originalstagemap,modifiedstagemap)
%A function that loads and compares two stage maps and returns the 
%linear number of stage positions which are different
%By Daniel Berger for MIT-BCS Seung / Harvard Lichtman, April 2010

txt=sprintf('Loading ORIGINAL stage map %s ...',originalstagemap); disp(txt);
[tree1, rootname1, dom1]=xml_read(originalstagemap);
txt=sprintf('Loading MANUALLY EDITED stage map %s ...',modifiedstagemap); disp(txt);
[tree2, rootname2, dom2]=xml_read(modifiedstagemap);

if (tree1.NumPoints~=tree2.NumPoints)
  disp('WARNING: Number of stage positions is different!');
end;

numpoints=min(tree1.NumPoints,tree2.NumPoints);

xcoords=zeros(numpoints,3);
ycoords=zeros(numpoints,3);
stagerots=zeros(numpoints,3);
beamrots=zeros(numpoints,3);
beam_wd=zeros(numpoints,3);
beam_stigx=zeros(numpoints,2);
beam_stigy=zeros(numpoints,2);
detector_b=zeros(numpoints,2);
detector_c=zeros(numpoints,2);

% Read out Coordinates
for p=1:1:tree1.NumPoints
  tag=sprintf('Ref%i',p);
  command=sprintf('x=tree1.%s.Stage.X;',tag); eval(command);
  command=sprintf('y=tree1.%s.Stage.Y;',tag); eval(command);
  command=sprintf('stagerot=tree1.%s.Stage.Rot;',tag); eval(command);
  command=sprintf('beamrot=tree1.%s.Beam.ScanRot;',tag); eval(command);
  xcoords(p,1)=x; ycoords(p,1)=y; stagerots(p,1)=stagerot; beamrots(p,1)=beamrot;
  command=sprintf('beam_wd(p,1)=tree1.%s.Beam.WD;',tag); eval(command);
  command=sprintf('beam_stigx(p,1)=tree1.%s.Beam.StigX;',tag); eval(command);
  command=sprintf('beam_stigy(p,1)=tree1.%s.Beam.StigY;',tag); eval(command);
  command=sprintf('detector_b(p,1)=tree1.%s.Detector.B;',tag); eval(command);
  command=sprintf('detector_c(p,1)=tree1.%s.Detector.C;',tag); eval(command);

  command=sprintf('x=tree2.%s.Stage.X;',tag); eval(command);
  command=sprintf('y=tree2.%s.Stage.Y;',tag); eval(command);
  command=sprintf('stagerot=tree2.%s.Stage.Rot;',tag); eval(command);
  command=sprintf('beamrot=tree2.%s.Beam.ScanRot;',tag); eval(command);
  xcoords(p,2)=x; ycoords(p,2)=y; stagerots(p,2)=stagerot; beamrots(p,2)=beamrot;
  command=sprintf('beam_wd(p,2)=tree2.%s.Beam.WD;',tag); eval(command);
  command=sprintf('beam_stigx(p,2)=tree2.%s.Beam.StigX;',tag); eval(command);
  command=sprintf('beam_stigy(p,2)=tree2.%s.Beam.StigY;',tag); eval(command);
  command=sprintf('detector_b(p,2)=tree2.%s.Detector.B;',tag); eval(command);
  command=sprintf('detector_c(p,2)=tree2.%s.Detector.C;',tag); eval(command);
end;

diffpos=(xcoords(:,1)~=xcoords(:,2)) | (ycoords(:,1)~=ycoords(:,2));
nrupdated=sum(diffpos);
txt=sprintf('%d updated stage positions detected.',nrupdated); disp(txt);