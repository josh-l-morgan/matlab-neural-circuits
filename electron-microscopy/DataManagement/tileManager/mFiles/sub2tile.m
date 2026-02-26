function[tileName,ID] = sub2tile(w,s,r,c)

tileName = sprintf('Tile_r%d-c%d_w%03.0f_sec%03.0f', r, c, w, s);


IDstr = sprintf('1%03.0f%03.0f%03.0f%03.0f%',w,s,r,c);
ID = str2num(IDstr);







