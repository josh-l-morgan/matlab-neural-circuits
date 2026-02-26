function[ID,w,s,r,c] = tile2ID(tileName)

%tileName = sprintf('Tile_r%d-c%d_%s_sec%03.0f', r, c, waf, sec);

rs = regexp(tileName,'_r');
cs = regexp(tileName,'-c');
ws = regexp(tileName,'_w');
secs = regexp(tileName,'_sec');

w = str2num(tileName(ws+2:secs-1));
s = str2num(tileName(secs+4:end));
r = str2num(tileName(rs+2:cs-1));
c = str2num(tileName(cs+2:ws-1));

IDstr = sprintf('1%03.0f%03.0f%03.0f%03.0f%',w,s,r,c);
ID = str2num(IDstr);






