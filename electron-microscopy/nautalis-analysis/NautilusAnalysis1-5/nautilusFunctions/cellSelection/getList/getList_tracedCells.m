function[cellList tracedAxons tracedTCR tracedLIN] = getList_tracedCells()

%%Make list of all traced thalamocortical cells and the rgcs that innervate
%%them

tracedTCR = [   67 108   109   116   117   120   123   129   133   134   148   156   159   162 ... 
   163   169   170   201   203   205   206   207   210   213   215   216 ...
   217   218   232   267   601   907   110   268   273   237   224 919 259];

%%By eye
tracedTCR = [tracedTCR 903 909 110  130 131 135 163 224 ];
tracedTCR = unique(tracedTCR);

%%
tracedLIN = [125]

%% subjectively chosen to have substantial tracing performed
 tracedAxons = [ 1001        1006        1009        1012        1014        1021        1023        1025        1027        1028        1029        1030        1031 ...
        1032        1033        1034        1036        1037        1040        1041        1050        1051        1053        1054        1055        1056 ...
        1058        1059        1060        2003        2004        2006        2007        2008        2009        2011        2016        2024        2027 ...
        2028        2030        2032        2033        2034        2035        2038        2041        5003        5004        5005        5102        5104 ...
        5106        5108        6000        6001        6101        7001        7004        7016        7018        7019        7021        7022        9101  9102];
tracedAxons = unique(tracedAxons);
%%

cellList = [tracedTCR tracedAxons tracedLIN];