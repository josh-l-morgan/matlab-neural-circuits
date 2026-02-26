    fid = fopen(filename, 'r');
filename = fopen(fid);
        fclose(fid);
        
        [format, fmt_s] = imftype(filename)
        fmt_s = imformats(format);
        
        
        #include "tiffio.h"
#include "tiff.h"
...

TIFFOpen('foo.tif', 'r');

TIFFClose(tif);

the mex file appears to compile OK:
mex readtiffstack.c -I/opt/include -L/opt/lib/libtiff.a