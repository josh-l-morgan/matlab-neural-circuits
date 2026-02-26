function[tileName,w,s,r,c] = tile2ID(ID)

IDstr = num2str(ID);
w = str2num(IDstr(2:4));
s = str2num(IDstr(5:7));
r = str2num(IDstr(8:10));
c = str2num(IDstr(11:13));

tileName = sprintf('Tile_r%d-c%d_w%03.0f_sec%03.0f', r, c, w, s);







