function result = istif(fileName)
%ISTIF: True (1) if file has '.tif' extension, wrong (0) otherwise

result = strcmp(fileName(end-3:end),'.tif');
