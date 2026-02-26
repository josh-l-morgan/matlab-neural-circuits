function param=xml_readparamfromfile(filename)
%Reads the image parameters from a XML metadata file as provided by the Lichtman lab.
%By Daniel Berger, July 13th 2009

%Read XML file and generate struct
data=xml2struct(filename);

%Read slice name
referenceinfo=xml_getchildfromname(data,'ReferenceInfo');
v=xml_getchildfromname(referenceinfo,'Name');
param.slicename=v.children.data;
stage=xml_getchildfromname(referenceinfo,'Stage');
v=xml_getchildfromname(stage,'Rot');
param.rot=str2num(v.children.data);

%Read tile proportions
tileinfo=xml_getchildfromname(data,'TileInfo');
v=xml_getchildfromname(tileinfo,'TileWidth');
param.tilewidth=str2num(v.children.data);
v=xml_getchildfromname(tileinfo,'TileHeight');
param.tileheight=str2num(v.children.data);
v=xml_getchildfromname(tileinfo,'NumTilesX');
param.nrofcolumns=str2num(v.children.data);
v=xml_getchildfromname(tileinfo,'NumTilesY');
param.nrofrows=str2num(v.children.data);

%Read magnification
v=xml_getchildfromname(data,'Width');
param.imgwidth_mu=str2num(v.children.data);
v=xml_getchildfromname(data,'Height');
param.imgheight_mu=str2num(v.children.data);
v=xml_getchildfromname(data,'PixelSize');
param.pixelsize_nm=str2num(v.children.data);

param.xpos=zeros(param.nrofrows,param.nrofcolumns);
param.ypos=zeros(param.nrofrows,param.nrofcolumns);
%param.starttime=zeros(param.nrofrows,param.nrofcolumns);

tiles=xml_getchildfromname(data,'Tiles');
for row=1:1:param.nrofrows
  for column=1:1:param.nrofcolumns
    tile=xml_getchildfromnameandattrib(tiles,'Tile',{'col', 'row'},{num2str(column), num2str(row)});
    v=xml_getchildfromname(tile,'StartTime');
    param.starttime{row,column}=v.children.data;
    v=xml_getchildfromname(tile,'StageX');
    param.xpos(row,column)=str2num(v.children.data);
    v=xml_getchildfromname(tile,'StageY');
    param.ypos(row,column)=str2num(v.children.data);
  end;
end;
