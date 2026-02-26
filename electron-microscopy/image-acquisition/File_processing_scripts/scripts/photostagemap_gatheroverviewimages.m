function param=photostagemap_gatheroverviewimages(param)
%This function is used to copy overview images to a local folder
%By Daniel Berger for MIT-BCS Seung / Harvard Lichtman, March 2010

lparam=param.overviewimages;

for w=param.processwafers
  nrofstrips=size(param.generatestagemap.stripslicepos{w},1);
  for st=1:1:nrofstrips
    nrofsections=size(param.generatestagemap.stripslicepos{w}{st},1);
    for sec=1:1:nrofsections
      sourcedir=sprintf(lparam.sourcedirectory,w);
      sourcesubdir=sprintf(lparam.sourcesubdirectory,w,st,sec);
      sourcename=sprintf(lparam.filenametemplate,w,st,sec);
      if (param.generatestagemap.stripslicepos{w}{st}(sec,3)==1) %this is a bad slice, add _b to the name
        sourcesubdir=[sourcesubdir(1:end-1) '_b' sourcesubdir(end)];
        sourcename=[sourcename(1:end-4) '_b' sourcename(end-3:end)];
      end;
      src=[sourcedir sourcesubdir sourcename];
      dst=[lparam.targetdirectory sourcename];
      txt=sprintf('Copying %s to %s ...',src,dst); disp(txt);
      [status,message,messageid]=copyfile(src,dst);
      if status==0
        disp('  ERROR: Copying failed.');
      end;

      
      %Get image parameters from the XML file
      xmlsourcename=sprintf(lparam.xmlnametemplate,w,st,sec);
      if (param.generatestagemap.stripslicepos{w}{st}(sec,3)==1) %this is a bad slice, add _b to the name
        xmlsourcename=[xmlsourcename(1:end-4) '_b' xmlsourcename(end-3:end)];
      end;
      xmlname=[sourcedir sourcesubdir xmlsourcename];
      [tree, rootname,dom]=xml_read(xmlname);
      param.overviewimages.stagex(w,st,sec)=tree.ReferenceInfo.Stage.X;
      param.overviewimages.stagey(w,st,sec)=tree.ReferenceInfo.Stage.Y;
      param.overviewimages.stagerot(w,st,sec)=tree.ReferenceInfo.Stage.Rot;
      param.overviewimages.scanrot(w,st,sec)=tree.ReferenceInfo.Beam.ScanRot;
      param.overviewimages.pixelsize_nm(w,st,sec)=tree.PixelSize.CONTENT;
      param.overviewimages.dwelltime_ns(w,st,sec)=tree.DwellTime.CONTENT;
      param.overviewimages.width_um(w,st,sec)=tree.Width.CONTENT;
      param.overviewimages.height_um(w,st,sec)=tree.Height.CONTENT;
      param.overviewimages.width_px(w,st,sec)=tree.TileInfo.TileWidth;
      param.overviewimages.height_px(w,st,sec)=tree.TileInfo.TileHeight;
      param.overviewimages.fov_um(w,st,sec)=tree.TileInfo.FOV.CONTENT;
    end;
  end;
end;