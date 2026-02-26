classdef VASTControlClass < handle
  
  %%%% BY DANIEL BERGER FOR HARVARD-LICHTMAN
  %%%% VERSION: 5.14, May 08, 2024 - for VAST Lite 1.5.0
  
  properties
    jtcpobj;
    jtcphelperclasspath;
    isconnected=0;
    
    inres=[];
    nrinints=0;
    inintdata=[];
    nrinuints=0;
    inuintdata=[];
    nrindoubles=0;
    indoubledata=[];
    nrinchars=0;
    inchardata=[];
    nrintext=0;
    intextdata={};
    nrinuint64s=0;
    inuint64data=[];
    
    parseheaderok=0;
    parseheaderlen=0;
    lasterror=0;
    thisversionnr=5; %Version 5 is for VAST Lite 1.3, 1.4, 1.5
    thissubversionnr=14;
    
    %islittleendian;
    indata=[];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTANTS

    GETINFO = 1;
    GETNUMBEROFSEGMENTS = 2;
    GETSEGMENTDATA     = 3;
    GETSEGMENTNAME     = 4;
    SETANCHORPOINT     = 5;
    SETSEGMENTNAME     = 6;
    SETSEGMENTCOLOR    = 7;
    GETVIEWCOORDINATES = 8;
    GETVIEWZOOM        = 9;
    SETVIEWCOORDINATES = 10;
    SETVIEWZOOM        = 11;
    GETNROFLAYERS      = 12;
    GETLAYERINFO       = 13;
    GETALLSEGMENTDATA  = 14;
    GETALLSEGMENTNAMES = 15;
    SETSELECTEDSEGMENTNR = 16;
    GETSELECTEDSEGMENTNR = 17;
    SETSELECTEDLAYERNR = 18;
    GETSELECTEDLAYERNR = 19;
    GETSEGIMAGERAW     = 20;
    GETSEGIMAGERLE     = 21;
    GETSEGIMAGESURFRLE = 22;
    SETSEGTRANSLATION  = 23;
    GETSEGIMAGERAWIMMEDIATE = 24;
    GETSEGIMAGERLEIMMEDIATE = 25;
    GETEMIMAGERAW      = 30;
    GETEMIMAGERAWIMMEDIATE = 31;
    REFRESHLAYERREGION = 32;
    GETPIXELVALUE      = 33;
    GETSCREENSHOTIMAGERAW = 40;
    GETSCREENSHOTIMAGERLE = 41;
    SETSEGIMAGERAW     = 50;
    SETSEGIMAGERLE     = 51;
    SETSEGMENTBBOX     = 60;
    GETFIRSTSEGMENTNR  = 61;
    GETHARDWAREINFO    = 62;
    ADDSEGMENT         = 63;
    MOVESEGMENT        = 64;
    GETDRAWINGPROPERTIES = 65;
    SETDRAWINGPROPERTIES = 66;
    GETFILLINGPROPERTIES = 67;
    SETFILLINGPROPERTIES = 68;
    
    GETANNOLAYERNROFOBJECTS = 70;
    GETANNOLAYEROBJECTDATA  = 71;
    GETANNOLAYEROBJECTNAMES  = 72;
    ADDNEWANNOOBJECT = 73;
    MOVEANNOOBJECT = 74;
    REMOVEANNOOBJECT = 75;
    SETSELECTEDANNOOBJECTNR = 76;
    GETSELECTEDANNOOBJECTNR = 77;
    GETAONODEDATA = 78;
    GETAONODELABELS = 79;
    
    SETSELECTEDAONODEBYDFSNR = 80;
    SETSELECTEDAONODEBYCOORDS = 81;
    GETSELECTEDAONODENR = 82;
    ADDAONODE           = 83;
    MOVESELECTEDAONODE          = 84;
    REMOVESELECTEDAONODE        = 85;
    SWAPSELECTEDAONODECHILDREN  = 86;
    MAKESELECTEDAONODEROOT      = 87;
    
    SPLITSELECTEDSKELETON = 88;
    WELDSKELETONS = 89;
    
    GETANNOOBJECT = 90;
    SETANNOOBJECT = 91;
    ADDANNOOBJECT = 92;
    GETCLOSESTAONODEBYCOORDS = 93;
    GETAONODEPARAMS = 94;
    SETAONODEPARAMS = 95;
    
    GETAPIVERSION      = 100;
    GETAPILAYERSENABLED = 101;
    SETAPILAYERSENABLED = 102;
    GETSELECTEDAPILAYERNR = 103;
    SETSELECTEDAPILAYERNR = 104;
    
    GETCURRENTUISTATE  = 110;
    GETERRORPOPUPSENABLED = 112;
    SETERRORPOPUPSENABLED = 113;
    SETUIMODE          = 114;
    SHOWWINDOW         = 115;
    SET2DVIEWORIENTATION = 116;
    GETPOPUPSENABLED = 117;
    SETPOPUPSENABLED = 118;
    
    ADDNEWLAYER        = 120;
    LOADLAYER          = 121;
    SAVELAYER          = 122;
    REMOVELAYER        = 123;
    MOVELAYER          = 124;
    SETLAYERINFO       = 125;
    GETMIPMAPSCALEFACTORS = 126;
    
    EXECUTEFILL        = 131;
    EXECUTELIMITEDFILL = 132;
    EXECUTECANVASPAINTSTROKE = 133;
    EXECUTESTARTAUTOSKELETONIZATION=134;
    EXECUTESTOPAUTOSKELETONIZATION=135;
    EXECUTEISAUTOSKELETONIZATIONDONE=136;
    
    SETTOOLPARAMETERS  = 151;
  end
    
  methods
    
    function obj=VASTControlClass()
      %Check for endianness of this computer
      %obj.islittleendian=0;
      %x=1; y=typecast(x,'uint8'); %[0 0 0 0 0 0 240 63] on little-endian
      %if (y(8)==63) obj.islittleendian=1; end;
      %javaaddpath(pwd);
      mp=mfilename('fullpath');
      m=mfilename;
      mpp=mp(1:end-size(m,2));
      javaaddpath(mpp);
      obj.jtcphelperclasspath=mpp;
    end
    
    function res=connect(obj,host,port,timeout)
      obj.jtcpobj=jtcp('request',host,port,'timeout',timeout,'serialize',false);
      if isfield(obj.jtcpobj,'error')
        res=0;
        return;
      end
      res=1; %If the request fails, a Java Error will stop the script, so there is no error if this line is reached (original version)
    end
    
    function res=disconnect(obj)
      jtcp('close',obj.jtcpobj);
      res=1; %If the request fails, a Java Error will stop the script, so there is no error if this line is reached
    end
    
    function res=getlasterror(obj)
      res=obj.lasterror;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GETINFO
    
    function [info, res] = getinfo(obj)
      obj.sendmessage(obj.GETINFO,[]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        if ((obj.nrinuints==7)&&(obj.nrindoubles==3)&&(obj.nrinints==3))
          info.datasizex=obj.inuintdata(1);
          info.datasizey=obj.inuintdata(2);
          info.datasizez=obj.inuintdata(3);
          info.voxelsizex=obj.indoubledata(1);
          info.voxelsizey=obj.indoubledata(2);
          info.voxelsizez=obj.indoubledata(3);
          info.cubesizex=obj.inuintdata(4);
          info.cubesizey=obj.inuintdata(5);
          info.cubesizez=obj.inuintdata(6);
          info.currentviewx=obj.inintdata(1);
          info.currentviewy=obj.inintdata(2);
          info.currentviewz=obj.inintdata(3);
          info.nrofmiplevels=obj.inuintdata(7);
        else
          info = [];
          res=0;
          obj.lasterror=2; %unexpected data
        end
      else
        info=[];
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2: getnumberofsegments();
    function [nr, res] = getnumberofsegments(obj)
      obj.sendmessage(obj.GETNUMBEROFSEGMENTS,[]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        if (obj.nrinuints==1)
          nr=obj.inuintdata(1);
        else
          nr = [];
          res=0;
          obj.lasterror=2; %unexpected data
        end
      else
        nr=0;
      end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3: //data=getsegmentdata(id);
    function [data, res] = getsegmentdata(obj, id)
      obj.sendmessage(obj.GETSEGMENTDATA,obj.bytesfromuint32(id));
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        if (obj.nrinints==9)&&(obj.nrinuints==9)
          %Columns: Nr  flags  red1 green1 blue1 pattern1  red2 green2 blue2 pattern2  anchorx anchory anchorz  parentnr childnr prevnr nextnr   collapsednr   bboxx1 bboxy1 bboxz1 bboxx2 bboxy2 bboxz2   "name"
          data.id=obj.inuintdata(1); %typecast(obj.inintdata(1),'uint32');
          data.flags=obj.inuintdata(2); %typecast(obj.inintdata(2),'uint32');
          data.col1=obj.inuintdata(3); %typecast(obj.inintdata(3),'uint32');
          data.col2=obj.inuintdata(4); %typecast(obj.inintdata(4),'uint32');
          data.anchorpoint=obj.inintdata(1:3);
          data.hierarchy=obj.inuintdata(5:8); %typecast(obj.inintdata(8:11),'uint32');
          data.collapsednr=obj.inuintdata(9); %typecast(obj.inintdata(12),'uint32');
          data.boundingbox=obj.inintdata(4:9);
        else
          data = [];
          res=0;
          obj.lasterror=2; %unexpected data received
        end
      else
        data=[];
      end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 4: getsegmentname(id)
    function [name, res] = getsegmentname(obj, id)

      obj.sendmessage(obj.GETSEGMENTNAME,obj.bytesfromuint32(id));
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        if (obj.inres==1)
          name = char(obj.inchardata);
        else
          name = [];
          res=0;
          obj.lasterror=2; %unexpected data received
        end
      else
        name=[];
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 5: setanchorpoint(id, x, y, z)
    function res = setanchorpoint(obj, id, x, y, z)
      obj.sendmessage(obj.SETANCHORPOINT,obj.bytesfromuint32([id x y z]));
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 6: setsegmentname(id, name)
    function res = setsegmentname(obj, id, name)
      obj.sendmessage(obj.SETSEGMENTNAME,[obj.bytesfromuint32(id) obj.bytesfromtext(name)]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 7: setsegmentcolor(id, r1,g1,b1,p1,r2,g2,b2,p2)
    function res = setsegmentcolor8(obj, id, r1,g1,b1,p1,r2,g2,b2,p2)
      %r1,g1,b1 is the primary color, r2,g2,b2 is the secondary color (all values 8bit 0..255). p1 is the pattern (0..16 allowed)
      
      v1=uint32(0);
      v2=v1;
      ir1=uint32(r1); ir2=uint32(r2);
      ig1=uint32(g1); ig2=uint32(g2);
      ib1=uint32(b1); ib2=uint32(b2);
      ip1=uint32(p1); ip2=uint32(p2);
      
      v1=v1+bitand(ip1,255)+bitshift(bitand(ib1,255),8)+bitshift(bitand(ig1,255),16)+bitshift(bitand(ir1,255),24);
      v2=v2+bitand(ip2,255)+bitshift(bitand(ib2,255),8)+bitshift(bitand(ig2,255),16)+bitshift(bitand(ir2,255),24);
      obj.sendmessage(obj.SETSEGMENTCOLOR,obj.bytesfromuint32([id v1 v2]));
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function res = setsegmentcolor32(obj, id, col1, col2)
      %col1 is the primary color, col2 is the secondary color (all values 32 bit).
      
      obj.sendmessage(obj.SETSEGMENTCOLOR,obj.bytesfromuint32([id uint32(col1) uint32(col2)]));
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GETVIEWCOORDINATES = 8;
    function [x,y,z, res] = getviewcoordinates(obj)
      %all coordinates in pixels at mip0

      obj.sendmessage(obj.GETVIEWCOORDINATES,[]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        if (obj.nrinints==3)
          x=obj.inintdata(1);
          y=obj.inintdata(2);
          z=obj.inintdata(3);        
        else
          x=0; y=0; z=0;
          res=0;
          obj.lasterror=2; %unexpected data received
        end
      else
        x=0; y=0; z=0;
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GETVIEWZOOM        = 9;
    function [zoom, res] = getviewzoom(obj)
      %all coordinates in pixels at mip0

      obj.sendmessage(obj.GETVIEWZOOM,[]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        if (obj.nrinints==1)
          zoom=obj.inintdata(1);
        else
          zoom=0;
          res=0;
          obj.lasterror=2; %unexpected data received
        end
      else
        zoom=0;
      end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SETVIEWCOORDINATES = 10;
    function res = setviewcoordinates(obj, x, y, z)
      %all coordinates in pixels at mip0

      obj.sendmessage(obj.SETVIEWCOORDINATES,obj.bytesfromuint32([uint32(x) uint32(y) uint32(z)]));
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SETVIEWZOOM        = 11;
    function res = setviewzoom(obj, zoom)
    %all coordinates in pixels at mip0

      obj.sendmessage(obj.SETVIEWZOOM,obj.bytesfromint32(int32(zoom)));
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GETNROFLAYERS = 12;
    function [nroflayers, res] = getnroflayers(obj)

      obj.sendmessage(obj.GETNROFLAYERS,[]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        if (obj.nrinuints==1)
          nroflayers=obj.inuintdata(1);
        else
          nroflayers=0;
          res=0;
          obj.lasterror=2; %unexpected data received
        end
      else
        nroflayers=0;
      end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GETLAYERINFO = 13;
    function [layerinfo, res] = getlayerinfo(obj, layernr)

      obj.sendmessage(obj.GETLAYERINFO,obj.bytesfromuint32(uint32(layernr)));
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        %if (obj.nrinints==8)&&((obj.nrinuints==2)||(obj.nrinuints==5))&&(obj.nrindoubles==3)
        if (obj.nrinints==8)&&(obj.nrinuints==7)&&(obj.nrindoubles==3)
          layerinfo.type=obj.inintdata(1);
          layerinfo.editable=obj.inintdata(2);
          layerinfo.visible=obj.inintdata(3);
          layerinfo.brightness=obj.inintdata(4);
          layerinfo.contrast=obj.inintdata(5);
          layerinfo.opacitylevel=obj.indoubledata(1);
          layerinfo.brightnesslevel=obj.indoubledata(2);
          layerinfo.contrastlevel=obj.indoubledata(3);
          layerinfo.blendmode=obj.inintdata(6);
          layerinfo.blendoradd=obj.inintdata(7);
          layerinfo.tintcolor=obj.inuintdata(1);
          layerinfo.name=obj.inchardata;
          %if (obj.nrinuints==5)
            layerinfo.redtargetcolor=obj.inuintdata(2);
            layerinfo.greentargetcolor=obj.inuintdata(3);
            layerinfo.bluetargetcolor=obj.inuintdata(4);
            layerinfo.bytesperpixel=obj.inuintdata(5);
          %else
          %end;
          layerinfo.ischanged=obj.inintdata(8);
          layerinfo.inverted=obj.inuintdata(6);
          layerinfo.solomode=obj.inuintdata(7);
        else
          layerinfo=[];
          res=0;
          obj.lasterror=2; %unexpected data received
        end
      else
        layerinfo=[];
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GETALLSEGMENTDATA  = 14;
    function [segdata, res] = getallsegmentdata(obj)

      obj.sendmessage(obj.GETALLSEGMENTDATA,[]);
      obj.readdatablockwithhelper();
      indata1=obj.indata;
      parseheader(obj,indata1); %parse header to find out how much data will be sent
      expectedlength=obj.parseheaderlen+12;
      obj.indata=int8(zeros(1,expectedlength));
      writepos=size(indata1,2)+1;
      obj.indata(1:size(indata1,2))=indata1;
      while (writepos<expectedlength)
        indata2=int8(jtcp('read',obj.jtcpobj,'helperClassPath',obj.jtcphelperclasspath));
        obj.indata(writepos:writepos+size(indata2,2)-1)=indata2;
        writepos=writepos+size(indata2,2);
      end
      
      if (obj.inres==0)
        parse(obj,obj.indata);
      end
      res=processerror(obj);
      if (res==1)
        uid=typecast(obj.indata(17:end), 'uint32');
        id=typecast(obj.indata(17:end), 'int32');
        sp=2;
        for i=1:1:uid(1)
          segdata{i}.id=i-1;
          segdata{i}.flags=uid(sp);
          segdata{i}.col1=uid(sp+1);
          segdata{i}.col2=uid(sp+2);
          segdata{i}.anchorpoint=id(sp+3:sp+5); % x,y,z
          segdata{i}.hierarchy=uid(sp+6:sp+9); % parent,child,prev,next
          segdata{i}.collapsednr=uid(sp+10);
          segdata{i}.boundingbox=id(sp+11:sp+16); %x1,y1,z1,x2,y2,z2
          
          sp=sp+17;
        end
      else
        segdata=[];
      end
    end
    
    
    function [segdatamatrix, res] = getallsegmentdatamatrix(obj)
      
      obj.sendmessage(obj.GETALLSEGMENTDATA,[]);
      obj.readdatablockwithhelper();
      indata1=obj.indata;
      parseheader(obj,indata1); %parse header to find out how much data will be sent
      expectedlength=obj.parseheaderlen+12;
      obj.indata=int8(zeros(1,expectedlength));
      writepos=size(indata1,2)+1;
      obj.indata(1:size(indata1,2))=indata1;
      while (writepos<expectedlength)
        indata2=int8(jtcp('read',obj.jtcpobj,'helperClassPath',obj.jtcphelperclasspath));
        obj.indata(writepos:writepos+size(indata2,2)-1)=indata2;
        writepos=writepos+size(indata2,2);
      end
      
      if (obj.inres==0)
        parse(obj,obj.indata);
      end
      res=processerror(obj);
      if (res==1)
        uid=typecast(obj.indata(17:end), 'uint32');
        id=typecast(obj.indata(17:end), 'int32');
        %bid=typecast(obj.indata(17:end), 'uint8'); %THIS DOES NOT WORK ??? WRONG RESULTS!
        sp=2;
        
        segdatamatrix=zeros(uid(1),24);
        
        for i=1:1:uid(1)
          segdatamatrix(i,1)=i-1;
          segdatamatrix(i,2)=uid(sp);
          segdatamatrix(i,3)=bitand(bitshift(uid(sp+1),-24),255);
          segdatamatrix(i,4)=bitand(bitshift(uid(sp+1),-16),255);
          segdatamatrix(i,5)=bitand(bitshift(uid(sp+1),-8),255);
          segdatamatrix(i,6)=bitand(uid(sp+1),255);
          segdatamatrix(i,7)=bitand(bitshift(uid(sp+2),-24),255);
          segdatamatrix(i,8)=bitand(bitshift(uid(sp+2),-16),255);
          segdatamatrix(i,9)=bitand(bitshift(uid(sp+2),-8),255);
          segdatamatrix(i,10)=bitand(uid(sp+2),255);
          segdatamatrix(i,11:13)=id(sp+3:sp+5);
          segdatamatrix(i,14:17)=uid(sp+6:sp+9);
          segdatamatrix(i,18)=uid(sp+10);
          segdatamatrix(i,19:24)=id(sp+11:sp+16);
          
          sp=sp+17;
        end
        segdatamatrix=segdatamatrix(2:end,:);
      else
        segdatamatrix=[];
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GETALLSEGMENTNAMES = 15;
    function [segname, res] = getallsegmentnames(obj)
      
      obj.sendmessage(obj.GETALLSEGMENTNAMES,[]);
      obj.readdatablockwithhelper();
      indata1=obj.indata;
      parseheader(obj,indata1); %parse header to find out how much data will be sent
      expectedlength=obj.parseheaderlen+12;
      obj.indata=int8(zeros(1,expectedlength));
      writepos=size(indata1,2)+1;
      obj.indata(1:size(indata1,2))=indata1;
      while (writepos<expectedlength)
        indata2=int8(jtcp('read',obj.jtcpobj,'helperClassPath',obj.jtcphelperclasspath));
        obj.indata(writepos:writepos+size(indata2,2)-1)=indata2;
        writepos=writepos+size(indata2,2);
      end
      
      if (obj.inres==0)
        parse(obj,obj.indata);
      end
      res=processerror(obj);
      if (res==1)
        nrnames=typecast(obj.indata(17:20), 'uint32');
        
        sp=20;
        for i=1:1:nrnames
          sq=sp+1;
          while (obj.indata(sq)~=0)
            sq=sq+1;
          end
          
          bf=obj.indata(sp+1:sq-1);
          segname{i}=char(bf);
          sp=sq;
        end
      else
        segname=[];
      end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SETSELECTEDSEGMENTNR = 16;
    function res=setselectedsegmentnr(obj,segmentnr)
      obj.sendmessage(obj.SETSELECTEDSEGMENTNR,obj.bytesfromint32(segmentnr));
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
      
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GETSELECTEDSEGMENTNR = 17;
    function [selectedsegmentnr, res]=getselectedsegmentnr(obj)
      
      obj.sendmessage(obj.GETSELECTEDSEGMENTNR,[]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        if (obj.nrinuints==1)
          selectedsegmentnr=obj.inuintdata(1);
        else
          selectedsegmentnr=-1;
          res=0;
          obj.lasterror=2; %unexpected data received
        end
      else
        selectedsegmentnr=-1;
      end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SETSELECTEDLAYERNR = 18;
    function res=setselectedlayernr(obj,layernr)
      obj.sendmessage(obj.SETSELECTEDLAYERNR,obj.bytesfromint32(layernr));
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GETSELECTEDLAYERNR  = 19;
    function [selectedlayernr, selectedemlayernr, selectedsegmentlayernr, res]=getselectedlayernr(obj)
      selectedlayernr=-1;
      selectedemlayernr=-1;
      selectedsegmentlayernr=-1;
      obj.sendmessage(obj.GETSELECTEDLAYERNR,[]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        if (obj.nrinints==3)
          selectedlayernr=obj.inintdata(1);
          selectedemlayernr=obj.inintdata(2);
          selectedsegmentlayernr=obj.inintdata(3);
        else
          if (obj.nrinints==5)
            selectedlayernr=obj.inintdata(1);
            selectedemlayernr=obj.inintdata(2);
            selectedsegmentlayernr=obj.inintdata(4);
          else
            res=0;
            obj.lasterror=2; %unexpected data received
          end
        end
      end
    end
    
    function [selectedlayernr, selectedemlayernr, selectedannolayernr, selectedsegmentlayernr, selectedtoollayernr, res]=getselectedlayernr2(obj)
      selectedlayernr=-1;
      selectedemlayernr=-1;
      selectedannolayernr=-1;
      selectedsegmentlayernr=-1;
      selectedtoollayernr=-1;
      obj.sendmessage(obj.GETSELECTEDLAYERNR,[]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        if (obj.nrinints==3) %for backwards compatibility
          selectedlayernr=obj.inintdata(1);
          selectedemlayernr=obj.inintdata(2);
          selectedsegmentlayernr=obj.inintdata(3);
        else
          if (obj.nrinints==5)
            selectedlayernr=obj.inintdata(1);
            selectedemlayernr=obj.inintdata(2);
            selectedannolayernr=obj.inintdata(3);
            selectedsegmentlayernr=obj.inintdata(4);
            selectedtoollayernr=obj.inintdata(5);
          else
            res=0;
            obj.lasterror=2; %unexpected data received
          end
        end
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GETSEGIMAGERAW = 20;
    % GETSEGIMAGERAWIMMEDIATE = 24;
    function [segimage, res] = getsegimageraw(obj, miplevel,minx,maxx,miny,maxy,minz,maxz, flipflag, immediateflag, requestloadflag)
      if ~exist('immediateflag','var') immediateflag=0; end
      if ~exist('requestloadflag','var') requestloadflag=0; end
      if (immediateflag==0)
        obj.sendmessage(obj.GETSEGIMAGERAW,obj.bytesfromuint32([uint32(miplevel) uint32(minx) uint32(maxx) uint32(miny) uint32(maxy) uint32(minz) uint32(maxz)]));
      else
        obj.sendmessage(obj.GETSEGIMAGERAWIMMEDIATE,obj.bytesfromuint32([uint32(miplevel) uint32(minx) uint32(maxx) uint32(miny) uint32(maxy) uint32(minz) uint32(maxz) uint32(requestloadflag)]));
      end
      obj.readdatablockwithhelper();
      indata1=obj.indata;
      parseheader(obj,indata1); %parse header to find out how much data will be sent
      expectedlength=obj.parseheaderlen+12;
      obj.indata=int8(zeros(1,expectedlength));
      writepos=size(indata1,2)+1;
      obj.indata(1:size(indata1,2))=indata1;
      while (writepos<expectedlength)
        indata2=int8(jtcp('read',obj.jtcpobj,'helperClassPath',obj.jtcphelperclasspath));
        obj.indata(writepos:writepos+size(indata2,2)-1)=indata2;
        writepos=writepos+size(indata2,2);
      end
      
      if (obj.inres==0)
        parse(obj,obj.indata);
      end
      res=processerror(obj);
      if (res==1)
        segimage=typecast(obj.indata(17:end), 'uint16');
        if (exist('flipflag','var'))
          if (flipflag==1)
            segimage=reshape(segimage,[maxx-minx+1,maxy-miny+1,maxz-minz+1]);
            segimage=permute(segimage,[2 1 3]);
            segimage=segimage(:);
          end
        end
      else
        segimage=[];
      end
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GETSEGIMAGERLE = 21; GETSEGIMAGESURFRLE = 22; GETSEGIMAGERLEIMMEDIATE = 25;
    function [segimageRLE,res] = getsegimageRLE(obj, miplevel,minx,maxx,miny,maxy,minz,maxz,surfonlyflag,immediateflag,requestloadflag)
      if ~exist('immediateflag','var') immediateflag=0; end
      if ~exist('requestloadflag','var') requestloadflag=0; end
      if (immediateflag==0)
        if (surfonlyflag==0)
          obj.sendmessage(obj.GETSEGIMAGERLE,obj.bytesfromuint32([uint32(miplevel) uint32(minx) uint32(maxx) uint32(miny) uint32(maxy) uint32(minz) uint32(maxz)]));
        else
          obj.sendmessage(obj.GETSEGIMAGESURFRLE,obj.bytesfromuint32([uint32(miplevel) uint32(minx) uint32(maxx) uint32(miny) uint32(maxy) uint32(minz) uint32(maxz)]));
        end
      else
        obj.sendmessage(obj.GETSEGIMAGERLEIMMEDIATE,obj.bytesfromuint32([uint32(miplevel) uint32(minx) uint32(maxx) uint32(miny) uint32(maxy) uint32(minz) uint32(maxz) uint32(requestloadflag)]));
      end
      obj.readdatablockwithhelper();
      indata1=obj.indata;
      parseheader(obj,indata1); %parse header to find out how much data will be sent
      expectedlength=obj.parseheaderlen+12;
      obj.indata=int8(zeros(1,expectedlength));
      writepos=size(indata1,2)+1;
      obj.indata(1:size(indata1,2))=indata1;
      while (writepos<expectedlength)
        indata2=int8(jtcp('read',obj.jtcpobj,'helperClassPath',obj.jtcphelperclasspath));
        obj.indata(writepos:writepos+size(indata2,2)-1)=indata2;
        writepos=writepos+size(indata2,2);
      end
      
      if (obj.inres==0)
        parse(obj,obj.indata);
      end
      res=processerror(obj);
      if (res==1)
        segimageRLE=typecast(obj.indata(17:end), 'uint16');
      else
        segimageRLE=[];
      end
    end

    
    function [segimage,res] = getsegimageRLEdecoded(obj, miplevel,minx,maxx,miny,maxy,minz,maxz,surfonlyflag, flipflag,immediateflag,requestloadflag)
      if ~exist('immediateflag','var') immediateflag=0; end
      if ~exist('requestloadflag','var') requestloadflag=0; end
      [segimageRLE, res] = getsegimageRLE(obj, miplevel,minx,maxx,miny,maxy,minz,maxz,surfonlyflag,immediateflag,requestloadflag);
      if (res==0)
        segimage=[];
      else
        %Decode RLE
        segimage=uint16(zeros(maxx-minx+1,maxy-miny+1,maxz-minz+1));
        dp=1;
     
        for sp=1:2:size(segimageRLE,2)
          val=segimageRLE(sp);
          num=segimageRLE(sp+1);
          segimage(dp:dp+uint32(num)-1)=val;
          dp=dp+uint32(num);
        end
      end
      
      if (exist('flipflag','var'))
        if (flipflag==1)
          segimage=permute(segimage,[2 1 3]);
        end
      end
    end
    
    
    function [values,numbers,res] = getRLEcountunique(obj, miplevel,minx,maxx,miny,maxy,minz,maxz,surfonlyflag,immediateflag,requestloadflag)
      if ~exist('immediateflag','var') immediateflag=0; end
      if ~exist('requestloadflag','var') requestloadflag=0; end
      [segimageRLE,res] = getsegimageRLE(obj, miplevel,minx,maxx,miny,maxy,minz,maxz,surfonlyflag,immediateflag,requestloadflag);
      if (res==0)
        values=[];
        numbers=[];
      else
        maxsegval=max(segimageRLE(1:2:end));
        na=zeros(maxsegval+1,1);
        %Decode RLE
     
        for sp=1:2:size(segimageRLE,2)
          val=segimageRLE(sp);
          num=segimageRLE(sp+1);
          na(val+1)=na(val+1)+double(num);
        end
        
% SOMEHOW THIS DOESNT WORK
%         sva=segimageRLE(1:2:end);
%         sna=segimageRLE(2:2:end);
%         na(sva+1)=na(sva+1)+double(sna)';
        
        values=find(na>0);
        numbers=na(values);
        values=values-1;
      end
    end
    
    
    function [segimage,values,numbers,res] = getsegimageRLEdecodedcountunique(obj, miplevel,minx,maxx,miny,maxy,minz,maxz,surfonlyflag, flipflag,immediateflag,requestloadflag)
      if ~exist('immediateflag','var') immediateflag=0; end
      if ~exist('requestloadflag','var') requestloadflag=0; end
      [segimageRLE,res] = getsegimageRLE(obj, miplevel,minx,maxx,miny,maxy,minz,maxz,surfonlyflag,immediateflag,requestloadflag);
      if (res==0)
        segimage=[];
        values=[];
        numbers=[];
      else
        maxsegval=max(segimageRLE(1:2:end));
        na=zeros(maxsegval+1,1);
        %Decode RLE
        segimage=uint16(zeros(maxx-minx+1,maxy-miny+1,maxz-minz+1));
        dp=1;
     
        for sp=1:2:size(segimageRLE,2)
          val=segimageRLE(sp);
          num=segimageRLE(sp+1);
          segimage(dp:dp+uint32(num)-1)=val;
          na(val+1)=na(val+1)+double(num);
          dp=dp+uint32(num);
        end
        
        values=find(na>0);
        numbers=na(values);
        values=values-1;
        
        if (exist('flipflag','var'))
          if (flipflag==1)
            segimage=permute(segimage,[2 1 3]);
          end
        end
      end
    end
    

    function [segimage,values,numbers,bboxes,res] = getsegimageRLEdecodedbboxes(obj, miplevel,minx,maxx,miny,maxy,minz,maxz,surfonlyflag, flipflag,immediateflag,requestloadflag)
      if ~exist('immediateflag','var') immediateflag=0; end
      if ~exist('requestloadflag','var') requestloadflag=0; end
      %tic
      [segimageRLE,res] = getsegimageRLE(obj, miplevel,minx,maxx,miny,maxy,minz,maxz,surfonlyflag,immediateflag,requestloadflag);
      %disp('image transmitted in:'); toc
      %tic
      if (res==0)
        segimage=[];
        values=[];
        numbers=[];
        bboxes=[];
      else
        maxsegval=max(segimageRLE(1:2:end));
        na=int32(zeros(maxsegval+1,1));
        bboxes=zeros(maxsegval+1,6)-1;
        %Decode RLE
        segimage=uint16(zeros(maxx-minx+1,maxy-miny+1,maxz-minz+1));
        dp=int32(1);
        xs=int32(maxx-minx+1); 
        ys=int32(maxy-miny+1); 
        zs=int32(maxz-minz+1);
        x1=1; y1=1; z1=1;
        for sp=1:2:size(segimageRLE,2)
          val=segimageRLE(sp);
          num=int32(segimageRLE(sp+1));
          dp2=dp+int32(num)-1;
          segimage(dp:dp2)=val;
          na(val+1)=na(val+1)+num;
          
          %%%%%%%%%%%%%% Bounding box computations
          if ((x1+num-1)<=xs)
            xmin=x1;
            xmax=x1+num-1;
            ymin=y1; ymax=y1;
            zmin=z1; zmax=z1;
            
            x1=x1+num; %go to next start
          else
          
            % Uses idivide because the standard division in matlab has the wrong rounding behavior
            z1=idivide(dp-1,xs*ys)+1; %z1=((dp-1)/(xs*ys))+1;
            r=dp-((z1-1)*xs*ys);
            y1=idivide(r-1,xs)+1; %y1=((r-1)/xs)+1;
            x1=r-((y1-1)*xs);
            
            z2=idivide(dp2-1,xs*ys)+1; %z2=((dp2-1)/(xs*ys))+1;
            r=dp2-((z2-1)*xs*ys);
            y2=idivide(r-1,xs)+1; %y2=((r-1)/xs)+1;
            x2=r-((y2-1)*xs);
            
            xmin=min([x1 x2]); xmax=max([x1 x2]);
            ymin=min([y1 y2]); ymax=max([y1 y2]);
            zmin=min([z1 z2]); zmax=max([z1 z2]);
            if (zmax>zmin)
              %we must go over the plane corner, which extends the bbox to max in XY
              xmin=1; xmax=xs;
              ymin=1; ymax=ys;
            end
            if (ymax>ymin)
              %we must go over the plane edge, which extends the bbox to max in X
              xmin=1; xmax=xs;
            end
            
            x1=x2+1; %go to next start
            y1=y2;
            z1=z2;
          end
          
          if (bboxes(val+1,1)==-1)
            bboxes(val+1,:)=[xmin,ymin,zmin,xmax,ymax,zmax];
          else
            
            tbbox=bboxes(val+1,:);
            if (xmin<tbbox(1))
              tbbox(1)=xmin;
            end
            if (ymin<tbbox(2))
              tbbox(2)=ymin;
            end
            if (zmin<tbbox(3))
              tbbox(3)=zmin;
            end
            if (xmax>tbbox(4))
              tbbox(4)=xmax;
            end
            if (ymax>tbbox(5))
              tbbox(5)=ymax;
            end
            if (zmax>tbbox(6))
              tbbox(6)=zmax;
            end
            bboxes(val+1,:)=tbbox;
          end
          
          dp=dp2+1;
        end
        
        values=find(na>0);
        numbers=na(values);
        bboxes=bboxes(values,:);
        values=values-1;
        
        if (exist('flipflag','var'))
          if (flipflag==1)
            segimage=permute(segimage,[2 1 3]);
          end
        end
      end
      %disp('  decoded in:');
      %toc
    end
    
    
    function tbbox=expandboundingbox(obj,bbox1,bbox2)
      %This assumes the following order of coordinates: xmin,ymin,zmin,xmax,ymax,zmax

      tbbox=bbox1;
      if (bbox2(1)<bbox1(1))
        tbbox(1)=bbox2(1);
      end
      if (bbox2(2)<bbox1(2))
        tbbox(2)=bbox2(2);
      end
      if (bbox2(3)<bbox1(3))
        tbbox(3)=bbox2(3);
      end
      if (bbox2(4)>bbox1(4))
        tbbox(4)=bbox2(4);
      end
      if (bbox2(5)>bbox1(5))
        tbbox(5)=bbox2(5);
      end
      if (bbox2(6)>bbox1(6))
        tbbox(6)=bbox2(6);
      end
    end
    
    function tbboxes=expandboundingboxes(obj,bboxes1,bboxes2)
      %This assumes the following order of coordinates: xmin,ymin,zmin,xmax,ymax,zmax
      %processes many bounding boxes at once
      %uses -1 in bboxes1 as indicator that the initial bounding box is defined by bboxes2
      
      tbboxes=bboxes1;
      for i=1:1:size(bboxes1,1)
        if (bboxes1(i,1)==-1)
          tbboxes(i,:)=bboxes2(i,:);
        else
          if (bboxes2(i,1)<bboxes1(i,1))
            tbboxes(i,1)=bboxes2(i,1);
          end
          if (bboxes2(i,2)<bboxes1(i,2))
            tbboxes(i,2)=bboxes2(i,2);
          end
          if (bboxes2(i,3)<bboxes1(i,3))
            tbboxes(i,3)=bboxes2(i,3);
          end
          if (bboxes2(i,4)>bboxes1(i,4))
            tbboxes(i,4)=bboxes2(i,4);
          end
          if (bboxes2(i,5)>bboxes1(i,5))
            tbboxes(i,5)=bboxes2(i,5);
          end
          if (bboxes2(i,6)>bboxes1(i,6))
            tbboxes(i,6)=bboxes2(i,6);
          end
        end
      end
    end
      

    function res = setsegtranslation(obj, sourcearray, targetarray)
      %sets the segmentation translation for getsegimage functions.
      %sourcearray is an array of segment numbers, targetarray an array of destination segment numbers.
      %before the image is transmitted, all voxels with a value in sourcearray will be translated to the corresponding number
      %in targetarray. segment numbers which do not appear in sourcearray will be set to 0 (background).
      %call with empty arrays to remove segmentation translation.

      if (length(sourcearray(:))~=length(targetarray(:)))
        %sourcearray and targetarray must have the same length.
        res=0;
        obj.lasterror=50;
      else
        translate=uint32(zeros(1,2*length(sourcearray(:))));
        translate(1:2:end)=uint32(sourcearray(:));
        translate(2:2:end)=uint32(targetarray(:));
        obj.sendmessage(obj.SETSEGTRANSLATION,obj.bytesfromdata(translate));
        obj.readdatablock();
        parse(obj,obj.indata);
        res=processerror(obj);
      end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GETEMIMAGERAW = 30;
    % GETEMIMAGERAWIMMEDIATE = 31;
    function [emimage,res] = getemimageraw(obj, layernr,miplevel,minx,maxx,miny,maxy,minz,maxz,immediateflag,requestloadflag)
      if ~exist('immediateflag','var') immediateflag=0; end
      if ~exist('requestloadflag','var') requestloadflag=0; end
      if (immediateflag==0)
        obj.sendmessage(obj.GETEMIMAGERAW,obj.bytesfromuint32([uint32(layernr) uint32(miplevel) uint32(minx) uint32(maxx) uint32(miny) uint32(maxy) uint32(minz) uint32(maxz)]));
      else
        obj.sendmessage(obj.GETEMIMAGERAWIMMEDIATE,obj.bytesfromuint32([uint32(layernr) uint32(miplevel) uint32(minx) uint32(maxx) uint32(miny) uint32(maxy) uint32(minz) uint32(maxz) uint32(requestloadflag)]));
        %obj.sendmessage(obj.GETEMIMAGERAWIMMEDIATE,obj.bytesfromuint32([uint32(layernr) uint32(miplevel) uint32(minx) uint32(maxx) uint32(miny) uint32(maxy) uint32(minz) uint32(maxz)]));
      end
      obj.readdatablockwithhelper();
      indata1=obj.indata;
      parseheader(obj,indata1); %parse header to find out how much data will be sent
      expectedlength=obj.parseheaderlen+12;
      obj.indata=int8(zeros(1,expectedlength));
      writepos=size(indata1,2)+1;
      obj.indata(1:size(indata1,2))=indata1;
      while (writepos<expectedlength)
        indata2=int8(jtcp('read',obj.jtcpobj,'helperClassPath',obj.jtcphelperclasspath));
        obj.indata(writepos:writepos+size(indata2,2)-1)=indata2;
        writepos=writepos+size(indata2,2);
      end
      
      if (obj.inres==0)
        parse(obj,obj.indata);
      end
      res=processerror(obj);
      if (res==1)
        emimage=typecast(obj.indata(17:end), 'uint8');
      else
        emimage=[];    
      end
    end
    
    function [emimage,res] = getemimage(obj, layernr,miplevel,minx,maxx,miny,maxy,minz,maxz,immediateflag,requestloadflag)
      if ~exist('immediateflag','var') immediateflag = 0; end
      if ~exist('requestloadflag','var') requestloadflag=0; end
      [emimageraw,res] = getemimageraw(obj, layernr,miplevel,minx,maxx,miny,maxy,minz,maxz,immediateflag,requestloadflag);
      bytesppx=int32(round(double(size(emimageraw,2))/((maxx-minx+1)*(maxy-miny+1)*(maxz-minz+1))));
      if (res==1)
        switch bytesppx
          case 1
            %One byte per pixel
            if (minz==maxz)
              emimage=permute(reshape(emimageraw,int32(maxx-minx+1),int32(maxy-miny+1)),[2 1]);
            else
              emimage=permute(reshape(emimageraw,int32(maxx-minx+1),int32(maxy-miny+1),int32(maxz-minz+1)),[2 1 3]);
            end
          case 3
            %Three bytes per pixel
            if (minz==maxz)
              emimage=flipdim(permute(reshape(emimageraw,3,int32(maxx-minx+1),int32(maxy-miny+1)),[3 2 1]),3);
            else
              emimage=flipdim(permute(reshape(emimageraw,3,int32(maxx-minx+1),int32(maxy-miny+1),int32(maxz-minz+1)),[3 2 4 1]),4);
            end
          case 4
            emimageraw=typecast(emimageraw,'uint32');
            if (minz==maxz)
              emimage=permute(reshape(emimageraw,int32(maxx-minx+1),int32(maxy-miny+1)),[2 1]);
            else
              emimage=permute(reshape(emimageraw,int32(maxx-minx+1),int32(maxy-miny+1),int32(maxz-minz+1)),[2 1 3]);
            end
          case 8
            emimageraw=typecast(emimageraw,'uint64');
            if (minz==maxz)
              emimage=permute(reshape(emimageraw,int32(maxx-minx+1),int32(maxy-miny+1)),[2 1]);
            else
              emimage=permute(reshape(emimageraw,int32(maxx-minx+1),int32(maxy-miny+1),int32(maxz-minz+1)),[2 1 3]);
            end
        end
      else
        emimage=emimageraw;
      end
    end
       

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % REFRESHLAYERREGION      = 32;
    function res = refreshlayerregion(obj, layernr, minx,maxx,miny,maxy,minz,maxz)
      obj.sendmessage(obj.REFRESHLAYERREGION,obj.bytesfromuint32([layernr minx maxx miny maxy minz maxz]));
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GETPIXELVALUE      = 33;
    function [pixelvalue, res]=getpixelvaluefromfullrescoords(obj, layernr, miplevel, x,y,z)
      pixelvalue=uint64(0);
      obj.sendmessage(obj.GETPIXELVALUE,obj.bytesfromuint32([layernr miplevel x y z]));
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        if (obj.nrinuint64s==1)
          pixelvalue=obj.inuint64data(1);
        else
          pixelvalue=0;
          res=0;
          obj.lasterror=2; %unexpected data received
        end
      else
        pixelvalue=0;
      end
    end
    
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GETSCREENSHOTIMAGERAW = 40;
    function [screenshotimage,res] = getscreenshotimageraw(obj, miplevel,minx,maxx,miny,maxy,minz,maxz,collapseseg)

      obj.sendmessage(obj.GETSCREENSHOTIMAGERAW,obj.bytesfromuint32([uint32(miplevel) uint32(minx) uint32(maxx) uint32(miny) uint32(maxy) uint32(minz) uint32(maxz) uint32(collapseseg)]));
      obj.readdatablockwithhelper();
      indata1=obj.indata;
      parseheader(obj,indata1); %parse header to find out how much data will be sent
      expectedlength=obj.parseheaderlen+12;
      obj.indata=int8(zeros(1,expectedlength));
      writepos=size(indata1,2)+1;
      obj.indata(1:size(indata1,2))=indata1;
      while (writepos<expectedlength)
        indata2=int8(jtcp('read',obj.jtcpobj,'helperClassPath',obj.jtcphelperclasspath));
        obj.indata(writepos:writepos+size(indata2,2)-1)=indata2;
        writepos=writepos+size(indata2,2);
      end
      
      if (obj.inres==0)
        parse(obj,obj.indata);
      end
      res=processerror(obj);
      if (res==1)
        screenshotimage=typecast(obj.indata(17:end), 'uint8');
      else
        screenshotimage=[];    
      end
    end
    
    % GETSCREENSHOTIMAGERLE = 41;
    function [screenshotimage,res] = getscreenshotimageRLE(obj, miplevel,minx,maxx,miny,maxy,minz,maxz,collapseseg)

      obj.sendmessage(obj.GETSCREENSHOTIMAGERLE,obj.bytesfromuint32([uint32(miplevel) uint32(minx) uint32(maxx) uint32(miny) uint32(maxy) uint32(minz) uint32(maxz) uint32(collapseseg)]));
      obj.readdatablockwithhelper();
      indata1=obj.indata;
      parseheader(obj,indata1); %parse header to find out how much data will be sent
      expectedlength=obj.parseheaderlen+12;
      obj.indata=int8(zeros(1,expectedlength));
      writepos=size(indata1,2)+1;
      obj.indata(1:size(indata1,2))=indata1;
      while (writepos<expectedlength)
        indata2=int8(jtcp('read',obj.jtcpobj,'helperClassPath',obj.jtcphelperclasspath));
        obj.indata(writepos:writepos+size(indata2,2)-1)=indata2;
        writepos=writepos+size(indata2,2);
      end
      
      if (obj.inres==0)
        parse(obj,obj.indata);
      end
      res=processerror(obj);
      if (res==1)
        rledata=typecast(obj.indata(17:end), 'uint8');
        
        if (length(rledata)==(maxx-minx+1)*(maxy-miny+1)*(maxz-minz+1)*3)
          screenshotimage=rledata;
        else
          %Decode RLE
          screenshotimage=uint8(zeros(1,(maxx-minx+1)*(maxy-miny+1)*(maxz-minz+1)*3));
          dp=1;
          
          for sp=1:4:size(rledata,2)
            val1=rledata(sp);
            val2=rledata(sp+1);
            val3=rledata(sp+2);
            num=rledata(sp+3);
            screenshotimage(dp:3:dp+uint32(num)*3-3)=val1;
            screenshotimage(dp+1:3:dp+1+uint32(num)*3-3)=val2;
            screenshotimage(dp+2:3:dp+2+uint32(num)*3-3)=val3;
            dp=dp+uint32(num)*3;
          end
        end
        
      else
        screenshotimage=[];    
      end
    end
    
    
    function [screenshotimage,res] = getscreenshotimage(obj, miplevel,minx,maxx,miny,maxy,minz,maxz,collapseseg,varargin)
      if (nargin>9)
        userle=varargin{1};
      else
        userle=0;
      end
      if (userle)
        [screenshotimageraw,res] = getscreenshotimageRLE(obj, miplevel,minx,maxx,miny,maxy,minz,maxz,collapseseg);
      else
        [screenshotimageraw,res] = getscreenshotimageraw(obj, miplevel,minx,maxx,miny,maxy,minz,maxz,collapseseg);
      end
      if (res==1)
        if (minz==maxz)
          screenshotimage=flipdim(permute(reshape(screenshotimageraw,3,int32(maxx-minx+1),int32(maxy-miny+1)),[3 2 1]),3);
        else
          screenshotimage=flipdim(permute(reshape(screenshotimageraw,3,int32(maxx-minx+1),int32(maxy-miny+1),int32(maxz-minz+1)),[3 2 4 1]),4);
        end
      else
        screenshotimage=screenshotimageraw;
      end
    end
    
    
    function res = orderscreenshotimage(obj, miplevel,minx,maxx,miny,maxy,minz,maxz,collapseseg)
      %Orders VAST to send data and returns immediately. Must call pickupscreenshotimage as next API function!
      obj.sendmessage(obj.GETSCREENSHOTIMAGERAW,obj.bytesfromuint32([uint32(miplevel) uint32(minx) uint32(maxx) uint32(miny) uint32(maxy) uint32(minz) uint32(maxz) uint32(collapseseg)]));
      res=1;
    end
      
    function [screenshotimage,res] = pickupscreenshotimage(obj, miplevel,minx,maxx,miny,maxy,minz,maxz,collapseseg)
      %Call after orderscreenshotimage to receive data. 
      %Parameters sent during order have to match parameters of pickup function!
      %Do not call other VASTControlClass functions between order and pickup!
      obj.readdatablockwithhelper();
      indata1=obj.indata;
      parseheader(obj,indata1); %parse header to find out how much data will be sent
      expectedlength=obj.parseheaderlen+12;
      obj.indata=int8(zeros(1,expectedlength));
      writepos=size(indata1,2)+1;
      obj.indata(1:size(indata1,2))=indata1;
      while (writepos<expectedlength)
        indata2=int8(jtcp('read',obj.jtcpobj,'helperClassPath',obj.jtcphelperclasspath));
        obj.indata(writepos:writepos+size(indata2,2)-1)=indata2;
        writepos=writepos+size(indata2,2);
      end
      
      if (obj.inres==0)
        parse(obj,obj.indata);
      end
      res=processerror(obj);
      if (res==1)
        screenshotimageraw=typecast(obj.indata(17:end), 'uint8');
      else
        screenshotimageraw=[];    
      end

      %[screenshotimageraw,res] = getscreenshotimageraw(obj, miplevel,minx,maxx,miny,maxy,minz,maxz,collapseseg);
      if (res==1)
        if (minz==maxz)
          screenshotimage=flipdim(permute(reshape(screenshotimageraw,3,int32(maxx-minx+1),int32(maxy-miny+1)),[3 2 1]),3);
        else
          screenshotimage=flipdim(permute(reshape(screenshotimageraw,3,int32(maxx-minx+1),int32(maxy-miny+1),int32(maxz-minz+1)),[3 2 4 1]),4);
        end
      else
        screenshotimage=screenshotimageraw;
      end
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SETSEGIMAGERAW = 50;
    function res = setsegimageraw(obj, miplevel,minx,maxx,miny,maxy,minz,maxz,segimage)

      if ((size(segimage,2)~=int32(maxx-minx+1))||(size(segimage,1)~=int32(maxy-miny+1))||(size(segimage,3)~=int32(maxz-minz+1)))
        %sourcearray and targetarray must have the same length.
        res=0;
        obj.lasterror=13;
      else
        mparams=obj.bytesfromuint32([uint32(miplevel) uint32(minx) uint32(maxx) uint32(miny) uint32(maxy) uint32(minz) uint32(maxz)]);
        segimage=permute(segimage,[2 1 3]);
        mdata=obj.bytesfromdata(uint16(segimage(:)));
        obj.sendmessage(obj.SETSEGIMAGERAW,[mparams mdata]);
        obj.readdatablock();
        parse(obj,obj.indata);
        res=processerror(obj);
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SETSEGIMAGERLE = 51;
    function res = setsegimageRLE(obj, miplevel,minx,maxx,miny,maxy,minz,maxz,segimage)
      %Encodes the image data as RLE and sends the encoded data to VAST
      if ((size(segimage,2)~=int32(maxx-minx+1))||(size(segimage,1)~=int32(maxy-miny+1))||(size(segimage,3)~=int32(maxz-minz+1)))
        %sourcearray and targetarray must have the same length.
        res=0;
        obj.lasterror=13;
      else
        mparams=obj.bytesfromuint32([uint32(miplevel) uint32(minx) uint32(maxx) uint32(miny) uint32(maxy) uint32(minz) uint32(maxz)]);
        segimage=permute(segimage,[2 1 3]);
        rledata=zeros(1,length(segimage),'uint16');

        % RLE encode
        rdp=1;
        rsp=1;
        val=segimage(rsp);
        num=1;
        rsp=rsp+1;
        dp=length(segimage(:));
        
        while ((rsp<=dp)&&(rdp<dp-2))
          if ((segimage(rsp)==val)&&(num<65535))
            num=num+1;
          else
            %store val,num pair here
            rledata(rdp)=val;
            rledata(rdp+1)=num;
            rdp=rdp+2;
            
            val=segimage(rsp);
            num=1;
          end
          rsp=rsp+1;
        end
        if (rsp==dp+1)
          rledata(rdp)=val;
          rledata(rdp+1)=num;
          rdp=rdp+2;
        end
        %%%%%%
        if (rsp<dp+1)
          % if RLE-encoded message is larger than original, fallback to RAW
          segimage=permute(segimage,[2 1 3]);
          res = setsegimageraw(obj, miplevel,minx,maxx,miny,maxy,minz,maxz,segimage);
        else
          mdata=obj.bytesfromdata(rledata(1:rdp-1));
          obj.sendmessage(obj.SETSEGIMAGERLE,[mparams mdata]);
          obj.readdatablock();
          parse(obj,obj.indata);
          res=processerror(obj);
        end
      end
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SETSEGMENTBBOX = 60;
    function res = setsegmentbbox(obj, id, minx,maxx,miny,maxy,minz,maxz)
      obj.sendmessage(obj.SETSEGMENTBBOX,obj.bytesfromuint32([id minx maxx miny maxy minz maxz]));
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GETFIRSTSEGMENTNR  = 61;
    function [firstsegmentnr, res]=getfirstsegmentnr(obj)
      obj.sendmessage(obj.GETFIRSTSEGMENTNR,[]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        if (obj.nrinuints==1)
          firstsegmentnr=obj.inuintdata(1);
        else
          firstsegmentnr=-1;
          res=0;
          obj.lasterror=2; %unexpected data received
        end
      else
        firstsegmentnr=-1;
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GETHARDWAREINFO    = 62;
    function [info, res] = gethardwareinfo(obj)
      obj.sendmessage(obj.GETHARDWAREINFO,[]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        if ((obj.nrinuints==1)&&(obj.nrindoubles==7)&&(obj.nrinints==0)&&(obj.nrintext==5))
          info.computername=obj.intextdata{1};
          info.processorname=obj.intextdata{2};
          info.processorspeed_ghz=obj.indoubledata(1);
          info.nrofprocessorcores=obj.inuintdata(1);
          info.tickspeedmhz=obj.indoubledata(2);
          info.mmxssecapabilities=obj.intextdata{3};
          info.totalmemorygb=obj.indoubledata(3);
          info.freememorygb=obj.indoubledata(4);
          info.graphicscardname=obj.intextdata{4};
          info.graphicsdedicatedvideomemgb=obj.indoubledata(5);
          info.graphicsdedicatedsysmemgb=obj.indoubledata(6);
          info.graphicssharedsysmemgb=obj.indoubledata(7);
          info.graphicsrasterizerused=obj.intextdata{5};
        else
          info = [];
          res=0;
          obj.lasterror=2; %unexpected data
        end
      else
        info=[];
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ADDSEGMENT = 63;
    function [id, res] = addsegment(obj, refid, nextorchild, name)
      %nextorchild: 0: next, 1: child
      obj.sendmessage(obj.ADDSEGMENT,[obj.bytesfromuint32([refid nextorchild]) obj.bytesfromtext(name)]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        if (obj.nrinuints==1)
          id=obj.inuintdata(1);
        else
          id=0;
          res=0;
          obj.lasterror=2; %unexpected data received
        end
      else
        id=0;
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MOVESEGMENT = 64;
    function res = movesegment(obj, id, refid, nextorchild)
      obj.sendmessage(obj.MOVESEGMENT,obj.bytesfromuint32([id refid nextorchild]));
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GETDRAWINGPROPERTIES = 65;
    function [props, res] = getdrawingproperties(obj)
      obj.sendmessage(obj.GETDRAWINGPROPERTIES,[]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        if ((obj.nrinuints==7)&&(obj.nrinints==1)&&(obj.nrindoubles==2))
          props.paintcursordiameter=obj.inuintdata(1);
          flags=obj.inuintdata(2);
          props.paintcursorlocked=bitand(flags,1);
          props.autofill=bitand(bitshift(flags,-1),1);
          props.zscrollenabled=bitand(bitshift(flags,-2),1);
          props.overwritemode=obj.inuintdata(3);
          props.mippaintrestriction=obj.inintdata(1);
          props.paintdepth=obj.inuintdata(4);
          flags=obj.inuintdata(5);
          props.useconditionalpainting=bitand(flags,1);
          props.cp_contiguousonly=bitand(bitshift(flags,-1),1);
          props.cp_method=obj.inuintdata(6);
          props.cp_sourcelayernr=obj.inuintdata(7);
          props.cp_lowvalue=obj.indoubledata(1);
          props.cp_highvalue=obj.indoubledata(2);
        else
          props=[];
          res=0;
          obj.lasterror=2; %unexpected data received
        end
      else
        props=[];
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SETDRAWINGPROPERTIES = 66;
    function res = setdrawingproperties(obj, props)
      xflags=uint32(0);
      msg=[];
      
      if (isfield(props,'paintcursorlocked'))
        msg=[msg obj.bytesfromint32(props.paintcursorlocked)];
        xflags=xflags+1;
      end
      if (isfield(props,'paintcursordiameter'))
        msg=[msg obj.bytesfromuint32(props.paintcursordiameter)];
        xflags=xflags+2;
      end
      if (isfield(props,'autofill'))
        msg=[msg obj.bytesfromint32(props.autofill)];
        xflags=xflags+4;
      end
      if (isfield(props,'zscrollenabled'))
        msg=[msg obj.bytesfromint32(props.zscrollenabled)];
        xflags=xflags+8;
      end
      if (isfield(props,'overwritemode'))
        msg=[msg obj.bytesfromuint32(props.overwritemode)];
        xflags=xflags+16;
      end
      if (isfield(props,'mippaintrestriction'))
        msg=[msg obj.bytesfromint32(props.mippaintrestriction)];
        xflags=xflags+32;
      end
      if (isfield(props,'paintdepth'))
        msg=[msg obj.bytesfromuint32(props.paintdepth)];
        xflags=xflags+64;
      end
      if (isfield(props,'useconditionalpainting'))
        msg=[msg obj.bytesfromint32(props.useconditionalpainting)];
        xflags=xflags+128;
      end
      if (isfield(props,'cp_contiguousonly'))
        msg=[msg obj.bytesfromint32(props.cp_contiguousonly)];
        xflags=xflags+256;
      end
      if (isfield(props,'cp_method'))
        msg=[msg obj.bytesfromuint32(props.cp_method)];
        xflags=xflags+512;
      end
      if (isfield(props,'cp_sourcelayernr'))
        msg=[msg obj.bytesfromint32(props.cp_sourcelayernr)];
        xflags=xflags+1024;
      end
      if (isfield(props,'cp_lowvalue'))
        msg=[msg obj.bytesfromdouble(props.cp_lowvalue)];
        xflags=xflags+2048;
      end
      if (isfield(props,'cp_highvalue'))
        msg=[msg obj.bytesfromdouble(props.cp_highvalue)];
        xflags=xflags+4096;
      end
      
      obj.sendmessage(obj.SETDRAWINGPROPERTIES,[obj.bytesfromuint32(xflags) msg]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GETFILLINGPROPERTIES = 67;
    function [props, res] = getfillingproperties(obj)
      obj.sendmessage(obj.GETFILLINGPROPERTIES,[]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        if ((obj.nrinuints==3)&&(obj.nrinints==2)&&(obj.nrindoubles==6))
          props.overwritemode=obj.inuintdata(1);
          props.mippaintrestriction=obj.inintdata(1);
          flags=obj.inuintdata(2);
          props.donotfillsourcecolorzero=bitand(flags,1);
          props.sourcelayersameastarget=bitand(bitshift(flags,-1),1);
          props.sourcelayernr=obj.inintdata(2);
          props.method=obj.inuintdata(3);
          props.lowvalue_x=obj.indoubledata(1);
          props.highvalue_x=obj.indoubledata(2);
          props.lowvalue_y=obj.indoubledata(3);
          props.highvalue_y=obj.indoubledata(4);
          props.lowvalue_z=obj.indoubledata(5);
          props.highvalue_z=obj.indoubledata(6);
        else
          props=[];
          res=0;
          obj.lasterror=2; %unexpected data received
        end
      else
        props=[];
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SETFILLINGPROPERTIES = 68;
    function res = setfillingproperties(obj, props)
      xflags=uint32(0);
      msg=[];
      
      if (isfield(props,'overwritemode'))
        msg=[msg obj.bytesfromint32(props.overwritemode)];
        xflags=xflags+1;
      end
      if (isfield(props,'mippaintrestriction'))
        msg=[msg obj.bytesfromint32(props.mippaintrestriction)];
        xflags=xflags+2;
      end
      if (isfield(props,'donotfillsourcecolorzero'))
        msg=[msg obj.bytesfromint32(props.donotfillsourcecolorzero)];
        xflags=xflags+4;
      end
      if (isfield(props,'sourcelayersameastarget'))
        msg=[msg obj.bytesfromint32(props.sourcelayersameastarget)];
        xflags=xflags+8;
      end
      if (isfield(props,'sourcelayernr'))
        msg=[msg obj.bytesfromint32(props.sourcelayernr)];
        xflags=xflags+16;
      end
      if (isfield(props,'method'))
        msg=[msg obj.bytesfromint32(props.method)];
        xflags=xflags+32;
      end

      if (isfield(props,'lowvalue_x'))
        msg=[msg obj.bytesfromdouble(props.lowvalue_x)];
        xflags=xflags+64;
      end
      if (isfield(props,'highvalue_x'))
        msg=[msg obj.bytesfromdouble(props.highvalue_x)];
        xflags=xflags+128;
      end
      if (isfield(props,'lowvalue_y'))
        msg=[msg obj.bytesfromdouble(props.lowvalue_y)];
        xflags=xflags+256;
      end
      if (isfield(props,'highvalue_y'))
        msg=[msg obj.bytesfromdouble(props.highvalue_y)];
        xflags=xflags+512;
      end
      if (isfield(props,'lowvalue_z'))
        msg=[msg obj.bytesfromdouble(props.lowvalue_z)];
        xflags=xflags+1024;
      end
      if (isfield(props,'highvalue_z'))
        msg=[msg obj.bytesfromdouble(props.highvalue_z)];
        xflags=xflags+2048;
      end
      
      obj.sendmessage(obj.SETFILLINGPROPERTIES,[obj.bytesfromuint32(xflags) msg]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GETANNOLAYERNROFOBJECTS = 70;
    function [nrofobjects, firstannoobjectnr, res] = getannolayernrofobjects(obj)
      obj.sendmessage(obj.GETANNOLAYERNROFOBJECTS,[]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        if (obj.nrinuints==2)
          nrofobjects=obj.inuintdata(1);
          firstannoobjectnr=obj.inuintdata(2);
        else
          nrofobjects=[];
          firstannoobjectnr=[];
          res=0;
          obj.lasterror=2; %unexpected data
        end
      else
        nrofobjects=0;
        firstannoobjectnr=[];
      end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GETANNOLAYEROBJECTDATA  = 71;
%     1 ID number of the annoobject
%     2 Type of the annoobject: 0: folder, 1: skeleton
%     3 Flags field of the segment as 16-bit value: currently bit 0 is 'isused', bit 8 is 'isexpanded'
%     4-7 Primary color as red, green, blue, pattern1
%     8-11 Secondary color as red, green, blue, pattern2 (unused)
%     12-14 XYZ coordinates of the annoobject's anchor point (in voxels)
%     15-18 IDs of parent, child, previous and next annoobject (0 if none)
%     19 If the annoobject is collapsed into a folder, this is the folder ID [currently unsupported]
%     20-25 Annoobject bounding box
%     26 Lower 32 bits of user flags
%     27 Higher 32 bits of user flags
%     28 Lower 32 bits of user ID
%     29 Higher 32 bits of user ID
    function [annolayerobjectdata, res] = getannolayerobjectdata(obj)
            
      obj.sendmessage(obj.GETANNOLAYEROBJECTDATA,[]);
      obj.readdatablockwithhelper();
      indata1=obj.indata;
      parseheader(obj,indata1); %parse header to find out how much data will be sent
      expectedlength=obj.parseheaderlen+12;
      obj.indata=int8(zeros(1,expectedlength));
      writepos=size(indata1,2)+1;
      obj.indata(1:size(indata1,2))=indata1;
      while (writepos<expectedlength)
        indata2=int8(jtcp('read',obj.jtcpobj,'helperClassPath',obj.jtcphelperclasspath));
        obj.indata(writepos:writepos+size(indata2,2)-1)=indata2;
        writepos=writepos+size(indata2,2);
      end
      
      if (obj.inres==0)
        parse(obj,obj.indata);
      end
      res=processerror(obj);
      if (res==1)
        uid=typecast(obj.indata(17:end), 'uint32');
        id=typecast(obj.indata(17:end), 'int32');
        sp=2;
        
        annolayerobjectdata=zeros(uid(1),29,'uint32');
        for i=1:1:uid(1)
          annolayerobjectdata(i,1)=i;
          annolayerobjectdata(i,2)=bitand(uid(sp),65535);
          annolayerobjectdata(i,3)=bitand(bitshift(uid(sp),-16),65535);
          annolayerobjectdata(i,4)=bitand(bitshift(uid(sp+1),-24),255);
          annolayerobjectdata(i,5)=bitand(bitshift(uid(sp+1),-16),255);
          annolayerobjectdata(i,6)=bitand(bitshift(uid(sp+1),-8),255);
          annolayerobjectdata(i,7)=bitand(uid(sp+1),255);
          annolayerobjectdata(i,8)=bitand(bitshift(uid(sp+2),-24),255);
          annolayerobjectdata(i,9)=bitand(bitshift(uid(sp+2),-16),255);
          annolayerobjectdata(i,10)=bitand(bitshift(uid(sp+2),-8),255);
          annolayerobjectdata(i,11)=bitand(uid(sp+2),255);
          annolayerobjectdata(i,12:14)=uid(sp+3:sp+5); %annoobject root coords are unsigned
          annolayerobjectdata(i,15:18)=uid(sp+6:sp+9); %hierarchy
          annolayerobjectdata(i,19)=uid(sp+10);
          annolayerobjectdata(i,20:25)=uid(sp+11:sp+16); 
          annolayerobjectdata(i,26:29)=uid(sp+17:sp+20); 
          %userflags are 64bit - not representable in a double. Use two 32-bit values
          %userids are 64bit - not representable in a double. Use two 32-bit values
          sp=sp+21;
        end
        %annolayerobjectdata=annolayerobjectdata(2:end,:);
      else
        annolayerobjectdata=[];
      end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GETANNOLAYEROBJECTNAMES  = 72;
    function [aoname, res] = getannolayerobjectnames(obj)
      aoname=[];
      obj.sendmessage(obj.GETANNOLAYEROBJECTNAMES,[]);
      obj.readdatablockwithhelper();
      indata1=obj.indata;
      parseheader(obj,indata1); %parse header to find out how much data will be sent
      expectedlength=obj.parseheaderlen+12;
      obj.indata=int8(zeros(1,expectedlength));
      writepos=size(indata1,2)+1;
      obj.indata(1:size(indata1,2))=indata1;
      while (writepos<expectedlength)
        indata2=int8(jtcp('read',obj.jtcpobj,'helperClassPath',obj.jtcphelperclasspath));
        obj.indata(writepos:writepos+size(indata2,2)-1)=indata2;
        writepos=writepos+size(indata2,2);
      end
      
      if (obj.inres==0)
        parse(obj,obj.indata);
      end
      res=processerror(obj);
      if (res==1)
        nrnames=typecast(obj.indata(17:20), 'uint32');

        sp=20;
        aoname=cell(nrnames,1);
        for i=1:1:nrnames
          sq=sp+1;
          while (obj.indata(sq)~=0)
            sq=sq+1;
          end
          
          bf=obj.indata(sp+1:sq-1);
          aoname{i}=char(bf);
          sp=sq;
        end
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ADDNEWANNOOBJECT = 73;
    function [id, res] = addnewannoobject(obj, refid, nextorchild, type, name)
      %type: 0: folder, 1: skeleton
      obj.sendmessage(obj.ADDNEWANNOOBJECT,[obj.bytesfromuint32([refid nextorchild type]) obj.bytesfromtext(name)]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        if (obj.nrinuints==1)
          id=obj.inuintdata(1);
        else
          id=0;
          res=0;
          obj.lasterror=2; %unexpected data received
        end
      else
        id=0;
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MOVEANNOOBJECT = 74;
    function res = moveannoobject(obj, id, refid, nextorchild)
      obj.sendmessage(obj.MOVEANNOOBJECT,obj.bytesfromuint32([id refid nextorchild]));
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % REMOVEANNOOBJECT = 75;
    function res = removeannoobject(obj, id)
      obj.sendmessage(obj.REMOVEANNOOBJECT,obj.bytesfromuint32([id]));
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SETSELECTEDANNOOBJECTNR = 76;
    function res = setselectedannoobjectnr(obj, id)
      obj.sendmessage(obj.SETSELECTEDANNOOBJECTNR,obj.bytesfromuint32([id]));
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GETSELECTEDANNOOBJECTNR = 77;
    function [id, res]=getselectedannoobjectnr(obj)
      obj.sendmessage(obj.GETSELECTEDANNOOBJECTNR,[]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        if (obj.nrinuints==1)
          id=obj.inuintdata(1);
        else
          id=-1;
          res=0;
          obj.lasterror=2; %unexpected data received
        end
      else
        id=-1;
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GETAONODEDATA = 78;
    function [aonodedata, res] = getaonodedata(obj)
      obj.sendmessage(obj.GETAONODEDATA,[]);
      obj.readdatablockwithhelper();
      indata1=obj.indata;
      parseheader(obj,indata1); %parse header to find out how much data will be sent
      expectedlength=obj.parseheaderlen+12;
      obj.indata=int8(zeros(1,expectedlength));
      writepos=size(indata1,2)+1;
      obj.indata(1:size(indata1,2))=indata1;
      while (writepos<expectedlength)
        indata2=int8(jtcp('read',obj.jtcpobj,'helperClassPath',obj.jtcphelperclasspath));
        obj.indata(writepos:writepos+size(indata2,2)-1)=indata2;
        writepos=writepos+size(indata2,2);
      end
      
      if (obj.inres==0)
        parse(obj,obj.indata);
      end
      res=processerror(obj);
      if (res==1)
        uid=typecast(obj.indata(17:end), 'uint32');
        id=typecast(obj.indata(17:end), 'int32');
        %did=typecast(obj.indata(17:end), 'double');
        sp=2;
        
        aonodedata=zeros(uid(1),14);
        for i=1:1:uid(1)
          aonodedata(i,1)=i-1;
          aonodedata(i,2)=bitand(uid(sp),255); %Flags, isselected
          aonodedata(i,3)=bitand(bitshift(uid(sp),-8),255); %Flags, edgeflags
          aonodedata(i,4)=bitand(bitshift(uid(sp),-16),255); %Flags, haslabel
          aonodedata(i,5)=bitand(bitshift(uid(sp),-24),255); %Flags, reserved
          vals=double(uid(sp+1:sp+6));
          vals(vals==4294967295)=-1; %vals(vals==hex2dec('FFFFFFFF'))=-1;
          aonodedata(i,6:11)=vals;
          idx=17+((i-1)*12+9)*4;
          aonodedata(i,12)=typecast(obj.indata(idx:idx+7),'double');
          aonodedata(i,13:14)=uid(sp+10:sp+11);

          sp=sp+12;
        end
        %annolayerobjectdata=annolayerobjectdata(2:end,:);
      else
        aonodedata=[];
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GETAONODELABELS = 79;
    function [nodenumbers, nodelabels, res] = getaonodelabels(obj)
            
      obj.sendmessage(obj.GETAONODELABELS,[]);
      obj.readdatablockwithhelper();
      indata1=obj.indata;
      parseheader(obj,indata1); %parse header to find out how much data will be sent
      expectedlength=obj.parseheaderlen+12;
      obj.indata=int8(zeros(1,expectedlength));
      writepos=size(indata1,2)+1;
      obj.indata(1:size(indata1,2))=indata1;
      while (writepos<expectedlength)
        indata2=int8(jtcp('read',obj.jtcpobj,'helperClassPath',obj.jtcphelperclasspath));
        obj.indata(writepos:writepos+size(indata2,2)-1)=indata2;
        writepos=writepos+size(indata2,2);
      end
      
      if (obj.inres==0)
        parse(obj,obj.indata);
      end
      res=processerror(obj);
      if (res==1)
        nroflabels=typecast(obj.indata(17:17+3), 'uint32');
        uid=typecast(obj.indata(21:21+nroflabels*4-1), 'uint32');

        nodenumbers=zeros(nroflabels,1);
        for i=1:nroflabels
          nodenumbers(i)=uid(i);
        end
        
        sp=16+(nroflabels+1)*4;
        for i=1:1:nroflabels
          sq=sp+1;
          while (obj.indata(sq)~=0)
            sq=sq+1;
          end
          
          bf=obj.indata(sp+1:sq-1);
          nodelabels{i}=char(bf);
          sp=sq;
        end
      else
        nodenumbers=[];
        nodelabels=[];
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SETSELECTEDAONODEBYDFSNR = 80;
    function res = setselectedaonodebydfsnr(obj, nr)
      obj.sendmessage(obj.SETSELECTEDAONODEBYDFSNR,obj.bytesfromuint32([nr]));
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SETSELECTEDAONODEBYCOORDS = 81;
    function res=setselectedaonodebycoords(obj,x,y,z)
      obj.sendmessage(obj.SETSELECTEDAONODEBYCOORDS,obj.bytesfromuint32([x y z]));
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GETSELECTEDAONODENR = 82;
    function [nr, res]=getselectedaonodenr(obj)
      obj.sendmessage(obj.GETSELECTEDAONODENR,[]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        if (obj.nrinuints==1)
          nr=obj.inuintdata(1);
        else
          nr=-1;
          res=0;
          obj.lasterror=2; %unexpected data received
        end
      else
        nr=-1;
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ADDAONODE           = 83;
    function res=addaonode(obj,x,y,z)
      obj.sendmessage(obj.ADDAONODE,obj.bytesfromuint32([x y z]));
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MOVESELECTEDAONODE          = 84;
    function res=moveselectedaonode(obj,x,y,z)
      obj.sendmessage(obj.MOVESELECTEDAONODE,obj.bytesfromuint32([x y z]));
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % REMOVESELECTEDAONODE        = 85;
    function res=removeselectedaonode(obj)
      obj.sendmessage(obj.REMOVESELECTEDAONODE,[]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SWAPSELECTEDAONODECHILDREN  = 86;
    function res=swapselectedaonodechildren(obj)
      obj.sendmessage(obj.SWAPSELECTEDAONODECHILDREN,[]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MAKESELECTEDAONODEROOT      = 87;
    function res=makeselectedaonoderoot(obj)
      obj.sendmessage(obj.MAKESELECTEDAONODEROOT,[]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SPLITSELECTEDSKELETON = 88;
    %Splits the selected skeleton by removing the parent edge of the node with depth-first-search number newrootdfsnr,
    %and generates a new object that contains the split-off tree part with that node being the root node of the new object
    function [id, res]=splitselectedskeleton(obj, newrootdfsnr, newname)
      obj.sendmessage(obj.SPLITSELECTEDSKELETON,[obj.bytesfromuint32([newrootdfsnr]) obj.bytesfromtext(newname)]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        if (obj.nrinuints==1)
          id=obj.inuintdata(1);
        else
          id=0;
          res=0;
          obj.lasterror=2; %unexpected data received
        end
      else
        id=0;
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % WELDSKELETONS = 89;
    % Welds tree of annoobjectnr2 to tree of annoobjectnr1 node tn with tree of node n.
    %//replaces tn with n, appends tn->aoparent to n->aoparent, and removes tn->aoparent
    function res=weldskeletons(obj, annoobjectnr1, nodedfsnr1, annoobjectnr2, nodedfsnr2)
      obj.sendmessage(obj.WELDSKELETONS,[obj.bytesfromuint32([annoobjectnr1 nodedfsnr1 annoobjectnr2 nodedfsnr2])]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
%       if (res==1)
%         if (obj.nrinuints==1)
%           id=obj.inuintdata(1);
%         else
%           id=0;
%           res=0;
%           obj.lasterror=2; %unexpected data received
%         end
%       else
%         id=0;
%       end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GETANNOOBJECT = 90;
    % Returns all the data of the annoobject with the given id number in the selected annotation layer.
    % If id is 0, the currently selected annoobject will be returned.
    function [aodata, res] = getannoobject(obj, id)
      obj.sendmessage(obj.GETANNOOBJECT,obj.bytesfromuint32(id));
      obj.readdatablockwithhelper();
      indata1=obj.indata;
      parseheader(obj,indata1); %parse header to find out how much data will be sent
      expectedlength=obj.parseheaderlen+12;
      obj.indata=int8(zeros(1,expectedlength));
      writepos=size(indata1,2)+1;
      obj.indata(1:size(indata1,2))=indata1;
      while (writepos<expectedlength)
        indata2=int8(jtcp('read',obj.jtcpobj,'helperClassPath',obj.jtcphelperclasspath));
        obj.indata(writepos:writepos+size(indata2,2)-1)=indata2;
        writepos=writepos+size(indata2,2);
      end
      
      if (obj.inres==0)
        parse(obj,obj.indata);
      end
      res=processerror(obj);
      
      if (res==1)
        [aodata, res]=decodetostruct(obj, obj.indata(17:end));
        if (res==0)
          aodata=[];
          obj.lasterror=2; %unexpected data received
        end
      else
        aodata=[];
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SETANNOOBJECT = 91;
    function [id, res] = setannoobject(obj, id, aodata)
      %aodata has to be a valid struct of the format as returned by getannoobject
      obj.sendmessage(obj.SETANNOOBJECT,[obj.bytesfromuint32(id) obj.bytesfromdata(encodefromstruct(obj, aodata))]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        if (obj.nrinuints==1)
          id=obj.inuintdata(1);
        else
          id=0;
          res=0;
          obj.lasterror=2; %unexpected data received
        end
      else
        id=0;
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ADDANNOOBJECT = 92;
    function [id, res] = addannoobject(obj, refid, nextorchild, aodata)
      %aodata has to be a valid struct of the format as returned by getannoobject
      obj.sendmessage(obj.ADDANNOOBJECT,[obj.bytesfromuint32([refid nextorchild]) obj.bytesfromdata(encodefromstruct(obj, aodata))]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        if (obj.nrinuints==1)
          id=obj.inuintdata(1);
        else
          id=0;
          res=0;
          obj.lasterror=2; %unexpected data received
        end
      else
        id=0;
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GETCLOSESTAONODEBYCOORDS = 93;
    function [annoobjectid, nodedfsnr, distance, res]=getclosestaonodebycoords(obj,x,y,z,maxdistance)
      obj.sendmessage(obj.GETCLOSESTAONODEBYCOORDS,obj.bytesfromuint32([x y z maxdistance]));
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      annoobjectid=-1;
      nodedfsnr=-1;
      distance=-1;
      if (res==1)
        if ((obj.nrinints==2)&&(obj.nrindoubles==1))
          annoobjectid=obj.inintdata(1);
          nodedfsnr=obj.inintdata(2);
          distance=obj.indoubledata(1);
        end
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GETAONODEPARAMS = 94;
    % Returns the data of the node at depth-first number nodedfsnr of the annoobject with the given id number in the selected annotation layer.
    % If id is 0, the currently selected annoobject will be returned.
    function [nodedata, res] = getaonodeparams(obj, annoobjectid, nodedfsnr)
      obj.sendmessage(obj.GETAONODEPARAMS,obj.bytesfromuint32([annoobjectid nodedfsnr]));
      obj.readdatablockwithhelper();
      indata1=obj.indata;
      parseheader(obj,indata1); %parse header to find out how much data will be sent
      expectedlength=obj.parseheaderlen+12;
      obj.indata=int8(zeros(1,expectedlength));
      writepos=size(indata1,2)+1;
      obj.indata(1:size(indata1,2))=indata1;
      while (writepos<expectedlength)
        indata2=int8(jtcp('read',obj.jtcpobj,'helperClassPath',obj.jtcphelperclasspath));
        obj.indata(writepos:writepos+size(indata2,2)-1)=indata2;
        writepos=writepos+size(indata2,2);
      end
      
      if (obj.inres==0)
        parse(obj,obj.indata);
      end
      res=processerror(obj);
      
      if (res==1)
        [nodedata, res]=decodetostruct(obj, obj.indata(17:end));
        if (res==0)
          nodedata=[];
          obj.lasterror=2; %unexpected data received
        end
      else
        nodedata=[];
      end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SETAONODEPARAMS = 95;
    function res = setaonodeparams(obj, annoobjectid, nodedfsnr, nodedata)
      %nodedata has to be a valid struct of the format as returned by getannonodeparams
      obj.sendmessage(obj.SETAONODEPARAMS,[obj.bytesfromuint32([annoobjectid nodedfsnr]) obj.bytesfromdata(encodefromstruct(obj, nodedata))]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
%       if (res==1)
%         if (obj.nrinuints==1)
%           id=obj.inuintdata(1);
%         else
%           id=0;
%           res=0;
%           obj.lasterror=2; %unexpected data received
%         end;
%       else
%         id=0;
%       end;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 100: getapiversion();
    function [version, res] = getapiversion(obj)
      obj.lasterror=0;
      obj.sendmessage(obj.GETAPIVERSION,[]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        if (obj.nrinuints==1)
          version=obj.inuintdata(1);
        else
          version = [];
          res=0;
          obj.lasterror=2; %unexpected data
        end
      else
        version=[];
      end
    end
    
    function [version, subversion, res] = getcontrolclassversion(obj)
      %Returns the locally defined version number of this script (VastControlClass.m)
      res=1;
      version=obj.thisversionnr;
      subversion=obj.thissubversionnr;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GETAPILAYERSENABLED = 101;
    function [isenabled, res] = getapilayersenabled(obj)
      obj.sendmessage(obj.GETAPILAYERSENABLED,[]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        %success
        if (obj.nrinuints==1)
          isenabled=obj.inuintdata(1);
        else
          isenabled=0;
        end
      else
        %command failed
        isenabled = [];
        res=0;
        obj.lasterror=2; %unexpected data
      end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SETAPILAYERSENABLED = 102;
    function res = setapilayersenabled(obj, enabledflag)
      obj.sendmessage(obj.SETAPILAYERSENABLED,obj.bytesfromuint32(enabledflag));
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GETSELECTEDAPILAYERNR = 103;
    function [selectedlayernr, selectedemlayernr, selectedannolayernr, selectedsegmentlayernr, selectedtoollayernr, res]=getselectedapilayernr(obj)
      selectedlayernr=-1;
      selectedemlayernr=-1;
      selectedannolayernr=-1;
      selectedsegmentlayernr=-1;
      selectedtoollayernr=-1;
      obj.sendmessage(obj.GETSELECTEDAPILAYERNR,[]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        if (obj.nrinints==5)
          selectedlayernr=obj.inintdata(1);
          selectedemlayernr=obj.inintdata(2);
          selectedannolayernr=obj.inintdata(3);
          selectedsegmentlayernr=obj.inintdata(4);
          selectedtoollayernr=obj.inintdata(5);
        else
          res=0;
          obj.lasterror=2; %unexpected data received
        end
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  SETSELECTEDAPILAYERNR = 104;
    function res=setselectedapilayernr(obj,layernr)
      obj.sendmessage(obj.SETSELECTEDAPILAYERNR,obj.bytesfromint32(layernr));
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  GETCURRENTUISTATE = 110;
    function [state, res] = getcurrentuistate(obj)
      obj.sendmessage(obj.GETCURRENTUISTATE,[]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        if ((obj.nrinuints==5)&&(obj.nrindoubles==0)&&(obj.nrinints==9)&&(obj.nrintext==0))
          state.mousecoordsx=obj.inintdata(1);
          state.mousecoordsy=obj.inintdata(2);
          state.lastleftclickx=obj.inintdata(3);
          state.lastleftclicky=obj.inintdata(4);
          state.lastleftreleasex=obj.inintdata(5);
          state.lastleftreleasey=obj.inintdata(6);
          state.mousecoordsz=obj.inintdata(7);
          state.clientwindowwidth=obj.inintdata(8);
          state.clientwindowheight=obj.inintdata(9);
          state.reservedflag=bitand(obj.inuintdata(1),1);
          state.lbuttondown=bitand(bitshift(obj.inuintdata(1),-1),1);
          state.rbuttondown=bitand(bitshift(obj.inuintdata(1),-2),1);
          state.mbuttondown=bitand(bitshift(obj.inuintdata(1),-3),1);
          state.ctrlpressed=bitand(bitshift(obj.inuintdata(1),-4),1);
          state.shiftpressed=bitand(bitshift(obj.inuintdata(1),-5),1);
          state.deletepressed=bitand(bitshift(obj.inuintdata(1),-6),1);
          state.spacepressed=bitand(bitshift(obj.inuintdata(1),-7),1);
          state.spacewaspressed=bitand(bitshift(obj.inuintdata(1),-8),1);
          state.uimode=obj.inuintdata(2);
          state.hoversegmentnr=obj.inuintdata(3);
          state.miplevel=obj.inuintdata(4);
          state.paintcursordiameter=obj.inuintdata(5);
        else
          info = [];
          res=0;
          obj.lasterror=2; %unexpected data
        end
      else
        info=[];
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  GETERRORPOPUPSENABLED = 112;
    function [isenabled, res] = geterrorpopupsenabled(obj, errornr)
      obj.sendmessage(obj.GETERRORPOPUPSENABLED,obj.bytesfromuint32(errornr));
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        %success
        if (obj.nrinuints==1)
          isenabled=obj.inuintdata(1);
        else
          isenabled=0;
        end
      else
        %command failed
        isenabled = [];
        res=0;
        obj.lasterror=2; %unexpected data
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  SETERRORPOPUPSENABLED = 113;
    function res = seterrorpopupsenabled(obj, errornr, enabledflag)
      obj.sendmessage(obj.SETERRORPOPUPSENABLED,obj.bytesfromuint32([errornr enabledflag]));
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  SETUIMODE = 114;
    function res = setuimode(obj, uimode)
      obj.sendmessage(obj.SETUIMODE,obj.bytesfromuint32([uimode]));
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SHOWWINDOW         = 115;
    function res = showwindow(obj, windownr, onoff, xpos, ypos, width, height)
      lonoff=1;
      lxpos=-1; lypos=-1;
      lw=0; lh=0;
      haspos=0; haswh=0;
      if (nargin>2) lonoff=onoff; end
      if (nargin>4) lxpos=xpos; lypos=ypos; haspos=1; end
      if (nargin>6) lw=width; lh=height; haswh=1; end
      obj.sendmessage(obj.SHOWWINDOW,obj.bytesfromint32([windownr, lonoff, haspos, lxpos, lypos, haswh, lw, lh]));
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SET2DVIEWORIENTATION = 116;
    function res = set2dvieworientation(obj, orientation, viewportnr)
      lorientation=0;
      lviewportnr=0;

      if (nargin>1) lorientation=orientation; end
      if (nargin>2) lviewportnr=viewportnr; end
      obj.sendmessage(obj.SET2DVIEWORIENTATION,obj.bytesfromint32([lorientation, lviewportnr]));
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  GETPOPUPSENABLED = 117;
    function [isenabled, res] = getpopupsenabled(obj, typenr)
      obj.sendmessage(obj.GETPOPUPSENABLED,obj.bytesfromuint32(typenr));
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        %success
        if (obj.nrinuints==1)
          isenabled=obj.inuintdata(1);
        else
          isenabled=0;
        end
      else
        %command failed
        isenabled = [];
        res=0;
        obj.lasterror=2; %unexpected data
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  SETPOPUPSENABLED = 118;
    function res = setpopupsenabled(obj, typenr, enabledflag)
      obj.sendmessage(obj.SETPOPUPSENABLED,obj.bytesfromuint32([typenr enabledflag]));
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  ADDNEWLAYER        = 120;
    function [id, res] = addnewlayer(obj, laytype, name, varargin)
      %laytype: 1: annotation, 2: segmentation (XY mipmapping), 3: tool, 4: segmentation with XYXYZ mipmapping (taken from EM layer)
      %varargin is an optional parameter. if specified, it defines refid (the layer nr after which this layer is to be inserted)
      %if not specified, the layer will be inserted after the selected layer or as the first and only layer
      if (nargin>3)
        refid=varargin{1};
      else
        refid=-1;
      end
      obj.sendmessage(obj.ADDNEWLAYER,[obj.bytesfromint32([laytype refid]) obj.bytesfromtext(name)]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        if (obj.nrinuints==1)
          id=obj.inuintdata(1);
        else
          id=0;
          res=0;
          obj.lasterror=2; %unexpected data received
        end
      else
        id=0;
      end
    end
    
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  LOADLAYER          = 121;
    function [id, res] = loadlayer(obj, filename, varargin)
      %varargin is an optional parameter. if specified, it defines refid (the layer nr after which this layer is to be inserted)
      %if not specified, the layer will be inserted after the selected layer or as the first and only layer
      if (nargin>2)
        refid=varargin{1};
      else
        refid=-1;
      end
      obj.sendmessage(obj.LOADLAYER,[obj.bytesfromint32(refid) obj.bytesfromtext(filename)]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        if (obj.nrinuints==1)
          id=obj.inuintdata(1);
        else
          id=0;
          res=0;
          obj.lasterror=2; %unexpected data received
        end
      else
        id=0;
      end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  SAVELAYER          = 122;
    function [res] = savelayer(obj, layernr, targetfilename, varargin)
      forceit=0;
      subformat=0;
      if (nargin>3)
        forceit=varargin{1};
      end
      if (nargin>4)
        subformat=varargin{2};
      end
      obj.sendmessage(obj.SAVELAYER,[obj.bytesfromuint32([layernr forceit subformat]) obj.bytesfromtext(targetfilename)]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  REMOVELAYER        = 123;
    function [res] = removelayer(obj, layernr, forceit)
      obj.sendmessage(obj.REMOVELAYER,[obj.bytesfromuint32([layernr forceit])]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  MOVELAYER          = 124;
    function [newlayernr, res] = movelayer(obj, movedlayernr, afterlayernr)
      obj.sendmessage(obj.MOVELAYER,[obj.bytesfromuint32([movedlayernr afterlayernr])]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        if (obj.nrinuints==1)
          id=obj.inuintdata(1);
        else
          id=0;
          res=0;
          obj.lasterror=2; %unexpected data received
        end
      else
        id=0;
      end
      newlayernr=id;
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  SETLAYERINFO = 125;
    function res = setlayerinfo(obj, layernr, layerinfo)
      xflags=uint32(0);
      msg=obj.bytesfromuint32(layernr);
      
      if (isfield(layerinfo,'editable'))
        msg=[msg obj.bytesfromint32(layerinfo.editable)];
        xflags=xflags+1;
      end
      if (isfield(layerinfo,'visible'))
        msg=[msg obj.bytesfromint32(layerinfo.visible)];
        xflags=xflags+2;
      end
      if (isfield(layerinfo,'brightness'))
        msg=[msg obj.bytesfromint32(layerinfo.brightness)];
        xflags=xflags+4;
      end
      if (isfield(layerinfo,'contrast'))
        msg=[msg obj.bytesfromint32(layerinfo.contrast)];
        xflags=xflags+8;
      end
      if (isfield(layerinfo,'opacitylevel'))
        msg=[msg obj.bytesfromdouble(layerinfo.opacitylevel)];
        xflags=xflags+16;
      end
      if (isfield(layerinfo,'brightnesslevel'))
        msg=[msg obj.bytesfromdouble(layerinfo.brightnesslevel)];
        xflags=xflags+32;
      end
      if (isfield(layerinfo,'contrastlevel'))
        msg=[msg obj.bytesfromdouble(layerinfo.contrastlevel)];
        xflags=xflags+64;
      end
      if (isfield(layerinfo,'blendmode'))
        msg=[msg obj.bytesfromint32(layerinfo.blendmode)];
        xflags=xflags+128;
      end
      if (isfield(layerinfo,'blendoradd'))
        msg=[msg obj.bytesfromint32(layerinfo.blendoradd)];
        xflags=xflags+256;
      end
      if (isfield(layerinfo,'tintcolor'))
        msg=[msg obj.bytesfromuint32(layerinfo.tintcolor)];
        xflags=xflags+512;
      end
      if (isfield(layerinfo,'redtargetcolor'))
        msg=[msg obj.bytesfromuint32(layerinfo.redtargetcolor)];
        xflags=xflags+1024;
      end
      if (isfield(layerinfo,'greentargetcolor'))
        msg=[msg obj.bytesfromuint32(layerinfo.greentargetcolor)];
        xflags=xflags+2048;
      end
      if (isfield(layerinfo,'bluetargetcolor'))
        msg=[msg obj.bytesfromuint32(layerinfo.bluetargetcolor)];
        xflags=xflags+4096;
      end
      if (isfield(layerinfo,'inverted'))
        msg=[msg obj.bytesfromuint32(layerinfo.inverted)];
        xflags=xflags+8192;
      end
      if (isfield(layerinfo,'solomode'))
        msg=[msg obj.bytesfromuint32(layerinfo.solomode)];
        xflags=xflags+16384;
      end
      
      obj.sendmessage(obj.SETLAYERINFO,[obj.bytesfromuint32(xflags) msg]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  GETMIPMAPSCALEFACTORS = 126;
    function [mipscalematrix, res] = getmipmapscalefactors(obj, layernr)
      obj.sendmessage(obj.GETMIPMAPSCALEFACTORS,[obj.bytesfromuint32(layernr)]);
      obj.readdatablockwithhelper();
      indata1=obj.indata;
      parseheader(obj,indata1); %parse header to find out how much data will be sent
      expectedlength=obj.parseheaderlen+12;
      obj.indata=int8(zeros(1,expectedlength));
      writepos=size(indata1,2)+1;
      obj.indata(1:size(indata1,2))=indata1;
      while (writepos<expectedlength)
        indata2=int8(jtcp('read',obj.jtcpobj,'helperClassPath',obj.jtcphelperclasspath));
        obj.indata(writepos:writepos+size(indata2,2)-1)=indata2;
        writepos=writepos+size(indata2,2);
      end
      
      if (obj.inres==0)
        parse(obj,obj.indata);
      end
      res=processerror(obj);
      if (res==1)
        uid=typecast(obj.indata(17:end), 'uint32');
        id=typecast(obj.indata(17:end), 'int32');
        sp=2;
        
        mipscalematrix=zeros(uid(1),3);
        
        for i=1:1:uid(1)
          mipscalematrix(i,1)=uid(sp);
          mipscalematrix(i,2)=uid(sp+1);
          mipscalematrix(i,3)=uid(sp+2);
          sp=sp+3;
        end
        mipscalematrix(1,:)=[];
      else
        mipscalematrix=[];
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  EXECUTEFILL        = 131;
    function res = executefill(obj, sourcelayernr, targetlayernr, x, y, z, mip)
      % x, y, z are seed coordinates in full resolution coordinates
      % mip is the mip level at which the fill is to be executed
      obj.sendmessage(obj.EXECUTEFILL,[obj.bytesfromuint32([sourcelayernr targetlayernr x y z mip])]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  EXECUTELIMITEDFILL = 132;
    function res = executelimitedfill(obj, sourcelayernr, targetlayernr, x, y, z, mip)
      % x, y, z are seed coordinates in full resolution coordinates
      % mip is the mip level at which the fill is to be executed
      obj.sendmessage(obj.EXECUTELIMITEDFILL,[obj.bytesfromuint32([sourcelayernr targetlayernr x y z mip])]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  EXECUTECANVASPAINTSTROKE = 133;
    function res = executecanvaspaintstroke(obj, coords)
      % coords is a matrix of n rows with 2 columns, representing a list of (x, y) mouse coordinates in window
      c2=coords';
      mdata=obj.bytesfromdata(uint32(c2(:)));
      obj.sendmessage(obj.EXECUTECANVASPAINTSTROKE,[obj.bytesfromuint32([size(coords,1) size(coords,2)]) mdata]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EXECUTESTARTAUTOSKELETONIZATION=134;
    function res = executestartautoskeletonization(obj, toollayernr, mip, nodedistance_mu, nodestep, regionpadding_mu)
      obj.sendmessage(obj.EXECUTESTARTAUTOSKELETONIZATION,[obj.bytesfromuint32([toollayernr mip nodestep]) obj.bytesfromdouble([nodedistance_mu regionpadding_mu])]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EXECUTESTOPAUTOSKELETONIZATION=135;
    function res = executestopautoskeletonization(obj)
      obj.sendmessage(obj.EXECUTESTOPAUTOSKELETONIZATION,[]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EXECUTEGETAUTOSKELETONIZATIONSTATE=136;
    function [state, res] = executegetautoskeletonizationstate(obj)
      obj.sendmessage(obj.EXECUTEISAUTOSKELETONIZATIONDONE,[]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      nrnodesadded=0;
      iswaiting=0;
      if (res==1)
        %success
        if (obj.nrinuints==3)
          isdone=obj.inuintdata(1);
          iswaiting=obj.inuintdata(2);
          nrnodesadded=obj.inuintdata(3);
        else
          isdone=0;
        end
      else
        %command failed
        isdone = [];
        res=0;
        obj.lasterror=2; %unexpected data
      end
      state.isdone=isdone;
      state.iswaiting=iswaiting;
      state.nrnodesadded=nrnodesadded;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SETTOOLPARAMETERS  = 151;
    function res = settoolparameters(obj, toollayernr, toolnodenr, params)
      %params has to be a valid struct of the format as returned by getannoobject
      obj.sendmessage(obj.SETTOOLPARAMETERS,[obj.bytesfromuint32([toollayernr toolnodenr]) obj.bytesfromdata(encodefromstruct(obj, params))]);
      obj.readdatablock();
      parse(obj,obj.indata);
      res=processerror(obj);
      if (res==1)
        if (obj.nrinuints==1)
          id=obj.inuintdata(1);
        else
          id=0;
          res=0;
          obj.lasterror=2; %unexpected data received
        end
      else
        id=0;
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % HELPER FUNCTIONS
    
    function ids=getidsfromexactname(obj, name)
      [namesarray, res]=obj.getallsegmentnames();
      if (res==0) 
        ids=[];
      else
        namesarray(1)=[];
        [dummy,ids]=ind2sub(size(namesarray), strmatch(name, namesarray, 'exact'));
      end
    end
    
    function ret=getimmediatechildids(obj,parentidlist)
      %Gets a list of the IDs of the segment's children's tree (if it exists)
      ret=[];
      [data, res]=obj.getallsegmentdatamatrix();

      if (res==1)
        pal=parentidlist(:);
        for p=1:1:size(pal,1)
          index=parentidlist(p);
          if (index==0)
            %Get all top-level IDs
            i=obj.getfirstsegmentnr();
            if (i>0)
              ret=[ret i];
              while (data(i,17)>0) %add sizes of all nexts
                i=data(i,17);
                ret=[ret i];
              end
            end
          else
            if (data(index,15)>0)
              i=data(index,15);
              ret=[ret i];
              while (data(i,17)>0) %add sizes of all nexts
                i=data(i,17);
                ret=[ret i];
              end
            end
          end
        end
      end
    end
    
    function ret=getchildtreeids2(obj,data,parentidlist)
      %Gets a list of the IDs of the segment's children's tree (if it exists)
      
      ret=[];
      %[data, res]=obj.getallsegmentdatamatrix();
      
      %if (res==1)
        pal=parentidlist(:);
        for p=1:1:size(pal,1)
          index=parentidlist(p);
          if (data(index,15)>0)
            i=data(index,15);
            ret=[ret i obj.getchildtreeids2(data,i)]; %Add size of child tree
            while (data(i,17)>0) %add sizes of all nexts
              i=data(i,17);
              ret=[ret i obj.getchildtreeids2(data,i)];
            end
          end
        end
      %end;
    end
    
    function ret=getchildtreeids(obj,parentidlist)
      %Gets a list of the IDs of the segment's children's tree (if it exists)
      
      ret=[];
      [data, res]=obj.getallsegmentdatamatrix();
      
      if (res==1)
        pal=parentidlist(:);
        for p=1:1:size(pal,1)
          index=parentidlist(p);
          if (data(index,15)>0)
            i=data(index,15);
            ret=[ret i obj.getchildtreeids2(data,i)]; %Add size of child tree
            while (data(i,17)>0) %add sizes of all nexts
              i=data(i,17);
              ret=[ret i obj.getchildtreeids2(data,i)];
            end
          end
        end
      end
    end
    
    function xyzsize=getdatasizeatmip(obj,mip,layernr)
      %Wrote this from scratch. I hope it's correct!
      xyzsize=[];
      info=obj.getinfo();
      msf=obj.getmipmapscalefactors(layernr);
      if ((mip>=0)&&(mip<info.nrofmiplevels)) %apparently info.nrofmiplevels counts mip0
        if (mip==0)
          xyzsize=double([info.datasizex info.datasizey info.datasizez]);
        else
          xyzsize=[floor(double(info.datasizex)/msf(mip,1)) floor(double(info.datasizey)/msf(mip,2)) floor(double(info.datasizez)/msf(mip,3))];
        end
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function res = processerror(obj)
      obj.lasterror=0;
      if (obj.inres==0)
        if (obj.nrinuints==1)
          obj.lasterror=obj.inuintdata(1); %get error number sent by Vast
        else
          obj.lasterror=1; %unknown error
        end
      end
      res=obj.inres;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function res=bytesfromint32(obj,value)
      val=int32(value);
      va=typecast(val,'uint8');
      res=[];
      for i=1:4:max(size(va))
        res=[res 4 va(i) va(i+1) va(i+2) va(i+3)];
      end
    end
    
    function res=bytesfromuint32(obj,value)
      val=uint32(value);
      va=typecast(val,'uint8');
      res=[];
      for i=1:4:max(size(va))
        res=[res 1 va(i) va(i+1) va(i+2) va(i+3)];
      end
    end
    
%     function res=bytesfromdouble(obj,value)
%       val=double(value(1,1));
%       res=typecast(val,'uint8');
% %       if (obj.islittleendian==0)
% %         ires=res;
% %         res=ires(end:-1:1);
% %       end;
%       res=[2 res];
%     end

    function res=bytesfromdouble(obj,value)
      val=double(value);
      va=typecast(val,'uint8');
      res=[];
      for i=1:8:max(size(va))
        res=[res 2 va(i) va(i+1) va(i+2) va(i+3) va(i+4) va(i+5) va(i+6) va(i+7)];
      end
    end

    
    function res=bytesfromtext(obj,value)
      val=uint8(value);
      res=[3 val 0];
    end
    
    function res=bytesfromdata(obj,value)
      val=typecast(value(:),'uint8');
      len=uint32(length(val));
      len=typecast(len,'uint8');
      res=[5 len(1) len(2) len(3) len(4) val'];
    end
    
    %--------------------------------------------------------------------
    function res=sendmessage(obj,messagenr,message)
      %The VAST Server expects data in the following format:
      %0..3: 'VAST', message header for binary messages
      %4..11: total size of data following; least-significant byte first
      %12..15: message number (uint32)
      %16...: message parameters (uint8!)
      
      len=uint64(max(size(message)))+4;
      len1=uint8(bitand(len,255));
      len2=uint8(bitand(bitshift(len,-8),255));
      len3=uint8(bitand(bitshift(len,-16),255));
      len4=uint8(bitand(bitshift(len,-24),255));
      len5=uint8(bitand(bitshift(len,-32),255));
      len6=uint8(bitand(bitshift(len,-40),255));
      len7=uint8(bitand(bitshift(len,-48),255));
      len8=uint8(bitand(bitshift(len,-56),255));
      mnr=uint32(messagenr);
      mnr1=uint8(bitand(mnr,255));
      mnr2=uint8(bitand(bitshift(mnr,-8),255));
      mnr3=uint8(bitand(bitshift(mnr,-16),255));
      mnr4=uint8(bitand(bitshift(mnr,-24),255));
      msg1=uint8(['VAST' len1 len2 len3 len4 len5 len6 len7 len8 mnr1 mnr2 mnr3 mnr4 message]);
      msg2=typecast(msg1, 'int8'); %necessary for lossless transmission of uint8 data
      jtcp('write',obj.jtcpobj,msg2);
    end
    
    %--------------------------------------------------------------------
    function res=sendmessagewithhelper(obj,messagenr,message)
      %The VAST Server expects data in the following format:
      %0..3: 'VAST', message header for binary messages
      %4..11: total size of data following; least-significant byte first
      %12..15: message number (uint32)
      %16...: message parameters (uint8!)
      
      len=uint64(max(size(message)))+4;
      len1=uint8(bitand(len,255));
      len2=uint8(bitand(bitshift(len,-8),255));
      len3=uint8(bitand(bitshift(len,-16),255));
      len4=uint8(bitand(bitshift(len,-24),255));
      len5=uint8(bitand(bitshift(len,-32),255));
      len6=uint8(bitand(bitshift(len,-40),255));
      len7=uint8(bitand(bitshift(len,-48),255));
      len8=uint8(bitand(bitshift(len,-56),255));
      mnr=uint32(messagenr);
      mnr1=uint8(bitand(mnr,255));
      mnr2=uint8(bitand(bitshift(mnr,-8),255));
      mnr3=uint8(bitand(bitshift(mnr,-16),255));
      mnr4=uint8(bitand(bitshift(mnr,-24),255));
      msg1=uint8(['VAST' len1 len2 len3 len4 len5 len6 len7 len8 mnr1 mnr2 mnr3 mnr4 message]);
      msg2=typecast(msg1, 'int8'); %necessary for lossless transmission of uint8 data
      jtcp('write',obj.jtcpobj,msg2,'helperClassPath',obj.jtcphelperclasspath);
    end
    
    %--------------------------------------------------------------------
    function readdatablock(obj)
      obj.indata=[];
      while (min(size(obj.indata))==0)
        obj.indata=jtcp('read',obj.jtcpobj);
        pause(0.001);
      end
    end
    
    %--------------------------------------------------------------------
    function readdatablockwithhelper(obj)
      obj.indata=[];
      while (min(size(obj.indata))==0)
        obj.indata=int8(jtcp('read',obj.jtcpobj,'helperClassPath',obj.jtcphelperclasspath));
        pause(0.001);
      end
    end
    
    %--------------------------------------------------------------------
    function parseheader(obj,indata)
      obj.parseheaderok=0;
      obj.parseheaderlen=0;
      obj.inres=[];
      if (max(size(indata))<16) 
        return;
      end
      header=typecast(indata(1:12), 'uint8');
      if ((header(1)~='V')||(header(2)~='A')||(header(3)~='S')||(header(4)~='T')) 
        return;
      end
      obj.parseheaderok=1;
      obj.parseheaderlen=(uint64(header(5)))+(bitshift(uint64(header(6)),8))+(bitshift(uint64(header(7)),16))+(bitshift(uint64(header(8)),24));
      obj.parseheaderlen=obj.parseheaderlen+(bitshift(uint64(header(9)),32))+(bitshift(uint64(header(10)),40))+(bitshift(uint64(header(11)),48))+(bitshift(uint64(header(12)),56));
      rf=indata(13:16);
      obj.inres = typecast(rf, 'int32');
    end
    
    %--------------------------------------------------------------------
    function parse(obj,indata)
      obj.nrinints=0;
      obj.inintdata=[];
      obj.nrinuints=0;
      obj.inuintdata=[];
      obj.nrindoubles=0;
      obj.indoubledata=[];
      obj.nrinchars=0;
      obj.inchardata=[];
      obj.nrintext=0;
      obj.intextdata={};
      obj.nrinuint64s=0;
      obj.inuint64data=[];
      obj.inres=[];
      if (max(size(indata))<12)
        return;
      end
      indata=typecast(indata, 'uint8');
      if ((indata(1)~='V')||(indata(2)~='A')||(indata(3)~='S')||(indata(4)~='T'))
        return;
      end
      len=(uint64(indata(5)))+(bitshift(uint64(indata(6)),8))+(bitshift(uint64(indata(7)),16))+(bitshift(uint64(indata(8)),24));
      len=len+(bitshift(uint64(indata(9)),32))+(bitshift(uint64(indata(10)),40))+(bitshift(uint64(indata(11)),48))+(bitshift(uint64(indata(12)),56));
      if (len~=(max(size(indata))-12))
        return;
      end
      
      rf=indata(13:16);
      obj.inres = typecast(rf, 'int32');      
      
      p=17;
      
      while (p<size(indata,2))
        switch indata(p)
          case 1 %unsigned int
            obj.nrinuints=obj.nrinuints+1;
            bf=indata(p+1:p+4);
            obj.inuintdata(obj.nrinuints) = typecast(bf, 'uint32');
            p=p+5;
            
          case 2 %double
            obj.nrindoubles=obj.nrindoubles+1;
            bf=indata(p+1:p+8);
            obj.indoubledata(obj.nrindoubles) = typecast(bf, 'double');
            p=p+9;
            
          case 3 %0-terminated text
            q=p+1;
            while ((q<size(indata,2))&&(indata(q)~=0))
              q=q+1;
            end
            if (indata(q)~=0)
              p=len+4;
            else
              bf=indata(p+1:q-1);
              obj.inchardata=char(bf);
              obj.nrinchars=max(size(bf));
              obj.nrintext=obj.nrintext+1;
              obj.intextdata{obj.nrintext}=char(bf);
              p=q+1;
            end
            
          case 4 %signed int
            obj.nrinints=obj.nrinints+1;
            bf=indata(p+1:p+4);
            obj.inintdata(obj.nrinints) = typecast(bf, 'int32');
            p=p+5;
            
          case 6 %uint64
            obj.nrinuint64s=obj.nrinuint64s+1;
            bf=indata(p+1:p+8);
            obj.inuint64data(obj.nrinuint64s) = typecast(bf, 'uint64');
            p=p+9;
            
          otherwise
            p=size(indata,2);
        end
      end
      obj.inintdata=int32(obj.inintdata);
      obj.inuintdata=uint32(obj.inuintdata);
    end
    
    %--------------------------------------------------------------------
    function [out, res]=decodetostruct(obj, indata)
      res=1;
      ptr=1;
      out=[];
      while (ptr<length(indata))
        %get variable name
        ptr2=ptr;
        while (indata(ptr2)~=0) ptr2=ptr2+1; end
        variablename=indata(ptr:ptr2-1);
        ptr=ptr2+1;
        %get variable type
        type=indata(ptr);
        ptr=ptr+1;
        switch type
          case 0
            %unsigned 32-bit integer value
            value=typecast(indata(ptr:ptr+3), 'uint32');
            ptr=ptr+4;
            eval(['out.' variablename ' = uint32(' sprintf('%u',value) ');']);
          case 1
            %signed 32-bit integer value
            value=typecast(indata(ptr:ptr+3), 'int32');
            ptr=ptr+4;
            eval(['out.' variablename ' = int32(' sprintf('%d',value) ');']);
          case 2
            %uint64 value
            value=typecast(indata(ptr:ptr+7), 'uint64');
            ptr=ptr+8;
            eval(['out.' variablename ' = ' sprintf('%u',value) ';']);
          case 3
            %double value
            value=typecast(indata(ptr:ptr+7), 'double');
            ptr=ptr+8;
            eval(['out.' variablename ' = ' sprintf('%f',value) ';']);
          case 4
            %character string
            ptr2=ptr;
            while (indata(ptr2)~=0) ptr2=ptr2+1; end
            value=indata(ptr:ptr2-1);
            ptr=ptr2+1;
            eval(['out.' variablename ' = ''' value ''';']);
          case 5
            %matrix of unsigned int values
            xsize=typecast(indata(ptr:ptr+3), 'uint32');
            ptr=ptr+4;
            ysize=typecast(indata(ptr:ptr+3), 'uint32');
            ptr=ptr+4;
            mtx=typecast(indata(ptr:ptr+xsize*ysize*4-1), 'uint32');
            mtx=reshape(mtx,xsize,ysize)';
            eval(['out.' variablename ' = mtx;']);
            ptr=ptr+xsize*ysize*4;
          case 6
            %matrix of int32 values
            xsize=typecast(indata(ptr:ptr+3), 'uint32');
            ptr=ptr+4;
            ysize=typecast(indata(ptr:ptr+3), 'uint32');
            ptr=ptr+4;
            mtx=typecast(indata(ptr:ptr+xsize*ysize*4-1), 'int32');
            mtx=reshape(mtx,xsize,ysize)';
            eval(['out.' variablename ' = mtx;']);
            ptr=ptr+xsize*ysize*4;
          case 7
            %matrix of uint64 values
            xsize=typecast(indata(ptr:ptr+3), 'uint32');
            ptr=ptr+4;
            ysize=typecast(indata(ptr:ptr+3), 'uint32');
            ptr=ptr+4;
            mtx=typecast(indata(ptr:ptr+xsize*ysize*8-1), 'uint64');
            mtx=reshape(mtx,xsize,ysize)';
            eval(['out.' variablename ' = mtx;']);
            ptr=ptr+xsize*ysize*8;
          case 8
            %matrix of double values
            xsize=typecast(indata(ptr:ptr+3), 'uint32');
            ptr=ptr+4;
            ysize=typecast(indata(ptr:ptr+3), 'uint32');
            ptr=ptr+4;
            mtx=typecast(indata(ptr:ptr+xsize*ysize*8-1), 'double');
            mtx=reshape(mtx,xsize,ysize)';
            eval(['out.' variablename ' = mtx;']);
            ptr=ptr+xsize*ysize*8;
          case 9
            %array of strings; make into cell array
            nrofstrings=typecast(indata(ptr:ptr+3), 'uint32');
            ptr=ptr+4;
            totaltextlength=typecast(indata(ptr:ptr+3), 'uint32');
            ptr=ptr+4;
            str=indata(ptr:ptr+totaltextlength-1);
            p=1;
            if (nrofstrings>0)
              for i=1:nrofstrings
                p2=p;
                while (str(p2)~=0) p2=p2+1; end
                ca{i}=char(str(p:p2-1));
                p=p2+1;
              end
              ptr=ptr+totaltextlength;
              eval(['out.' variablename ' = ca;']);
            end
        end
      end
    end
    
    function [out, res]=encodefromstruct(obj, indata)
      res=1;
      if (~isa(indata,'struct'))
        res=0;
        out=[];
      else
        fields=fieldnames(indata);
        out=[];
        for i=1:length(fields)
          value = getfield(indata,fields{i});
          if ((min(size(value))==1)&&(max(size(value))==1)&&(~isa(value,'cell')))
            %This is a single-value field
            if (isa(value,'uint32'))
              out=[out unicode2native(fields{i}) 0 0 typecast(value, 'uint8')];
            end
            if (isa(value,'int32'))
              out=[out unicode2native(fields{i}) 0 1 typecast(value, 'uint8')];
            end
            if (isa(value,'int64'))
              out=[out unicode2native(fields{i}) 0 2 typecast(value, 'uint8')];
            end
            if (isa(value,'double'))
              out=[out unicode2native(fields{i}) 0 3 typecast(value, 'uint8')];
            end
          else
            %This is an array or matrix
            if (isa(value,'char'))
              out=[out unicode2native(fields{i}) 0 4 unicode2native(value) 0];
            end
            if (isa(value,'uint32'))
              mtx=value';
              xsize=uint32(size(mtx,1));
              ysize=uint32(size(mtx,2));
              out=[out unicode2native(fields{i}) 0 5 typecast(xsize,'uint8') typecast(ysize,'uint8') typecast(mtx(:)', 'uint8')];
            end
            if (isa(value,'int32'))
              mtx=value';
              xsize=uint32(size(mtx,1));
              ysize=uint32(size(mtx,2));
              out=[out unicode2native(fields{i}) 0 6 typecast(xsize,'uint8') typecast(ysize,'uint8') typecast(mtx(:)', 'uint8')];
            end
            if (isa(value,'uint64'))
              mtx=value';
              xsize=uint32(size(mtx,1));
              ysize=uint32(size(mtx,2));
              out=[out unicode2native(fields{i}) 0 7 typecast(xsize,'uint8') typecast(ysize,'uint8') typecast(mtx(:)', 'uint8')];
            end
            if (isa(value,'double'))
              mtx=value';
              xsize=uint32(size(mtx,1));
              ysize=uint32(size(mtx,2));
              out=[out unicode2native(fields{i}) 0 8 typecast(xsize,'uint8') typecast(ysize,'uint8') typecast(mtx(:)', 'uint8')];
            end
            if (isa(value,'cell'))
              %assued to be an array of strings in a 1-dim cell array
              cstr=[];
              nrofstrings=uint32(length(value));
              for s=1:length(value)
                cstr=[cstr unicode2native(value{s}) 0];
              end
              totaltextlength=uint32(length(cstr));
              out=[out unicode2native(fields{i}) 0 9 typecast(nrofstrings,'uint8') typecast(totaltextlength,'uint8') cstr];
            end
          end
        end
      end
    end
  end
end