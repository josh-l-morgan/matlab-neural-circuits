%This is a meta-script, which runs other scripts in the processing pipeline 
%for EM image alignment
%By Daniel Berger, June 1st 2009

%This script is written for the hip29nm Hippocampus data.

%Specification of stack parameters. These parameters have to be specified by the user for each new stack
parameterfile='W2W7_parameters.mat'; %this .mat file stores all information for the alignment
finecorrparamfile='W2W7_finecorrparam.mat'; %extra parameter file for the huge results of the fine-grain crosscorrelation

if exist(parameterfile,'file')
  load(parameterfile);
end;

%param.origdir='../../original/';
param.rawdir='../../raw/';
param.scaleddir='../../scaled/';
%param.baseorigname='_W%03dS%03dR%dC%d.tif';
param.baserawname='ROIb_W2-W7_s%03d.tif';
param.basescaledname='ROIb_W2-W7_s%03d_scaled16.png';
param.nrofslices=305;
param.firstslice=1; %The number of the first slice in the stack
param.nrofrows=1;
param.nrofcolumns=1;

save(parameterfile,'param');

%These are flags for different parts of the processing pipeline
flag.unwafer                =0;   %To convert images in W*S*R*C* file naming to S*R*C*
flag.reorder                =0;   %Image reordering and inverting
flag.checkimages            =0;   %Check whether images are correct
flag.writescanfaultfile     =0;   %Writes out a text file which contains a report of scan faults
flag.downscale              =0;   %Downscale images and store
flag.absrot                 =0;   %Compute absolute rigid orientation from cutting marks
flag.correctabsrotbetween   =0;   %Find and correct absolute rotation estimate outliers by comparing across slices
flag.correctabsrotwithin    =0;   %Find and correct absolute rotation estimate outliers by comparing tiles within a slice
flag.abstransbetweenslices  =0;   %Compute relative rigid translations between slices, based on absrot
flag.abstransneighbors      =0;   %Compute relative rigid translations between tiles, based on absrot
flag.relrot                 =0;   %Compute relative rigid rotations/translations by image comparison (for same-tile between-slices)
flag.relrotcheckconsistency =0;
flag.refinerelrot           =0;   %Refine the relative rigid rotations/translations with higher resolution (for same-tile between-slices)
flag.reviewrigid            =0;   %View the results of the rigid alignment
flag.computerigidmatrices   =0;   %Compute relative and absolute rigid transformation matrices
flag.renderrigid            =0;   %Render out a stack of rigidly aligned images for validation

flag.predictcorr            =0;  %Predict where corresponding points should be on neighboring slices, based on rigid alignment
flag.computelocalcorr       =0;  %Compute actual peak of cross-correlation between local image regions
flag.fititerativeaffine     =0;  %Fit affine transformations to the corresponding pairs found by computelocalcorr
flag.renderaffine           =1;  %Render out a stack of affinely aligned images
flag.renderaffineregion     =0;  %Render out a cropped region of the stack

flag.computefinelocalcorr   =0;  %Compute actual peak of cross-correlation between local image regions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Images un-wafer from directory 'original' to 'raw'
if flag.unwafer==1
  disp('---------------------------------------------------------------------');
  disp('Un-Wafering images ...');
  disp('---------------------------------------------------------------------');
  param.unwafer.nrofwafers=6;
  param.unwafer.nrofslices=[58 56 55 47 49 48];
  param.unwafer.targetfileformat='png';
  param.unwafer.invert=1;
  param=pipeline_unwafer(param);
  save(parameterfile,'param'); %this also computes param.nrofslices new
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image reordering
if flag.reorder==1
  disp('---------------------------------------------------------------------');
  disp('Reordering images ...');
  disp('---------------------------------------------------------------------');
%   param.reorder.sourceslice=[8 7 6 5 4 3 2 1 19 18  17 16 15 14 13 12 11 10 9 29 ...
%                               28 27 26 25 24 23 22 21 20 39  38 37 36 35 34 33 32 31 30 49 ...
%                               48 47 46 45 44 43 42 41 40];
%  param.reorder.sourceslice=[1 2 3 4 5 6 7 8 9 10  11 12 13 14 15 16 17 18 19 20 ...
%                             21 22 23 24 25 26 27 28 29 30  31 32 33 34 35 36 37 38 39 40 ...
%                             41 42 43 44 45 46 47 48 49];
%   param.reorder.sourceslice=[10 9 8 7 6 5 4 3 2 1  19 18 17 16 15 14 13 12 11 29 ...
%                              28 27 26 25 24 23 22 21 20 39  38 37 36 35 34 33 32 31 30 48 ...
%                              47 46 45 44 43 42 41 40];
  param.reorder.invert=[0 0 0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0 0 0 ...
                        0 0 0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0 0 0 ... 
                        0 0 0 0 0 0 0 0 0]; %on target ordering, not source!
  param.reorder.targetfileformat='png';
  param.reorder.dopart=0;       %if this is 1, only part of the images are redone (from firstslice to lastslice)
  param.reorder.firstslice=37;
  param.reorder.lastslice=39;
  param.reorder.targetfirstslice=1;
  pipeline_invert_reorder(param);
  save(parameterfile,'param');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking images for scan faults
if flag.checkimages==1
  disp('---------------------------------------------------------------------');
  disp('Checking images for scan faults ...');
  disp('---------------------------------------------------------------------');
  param.checkimages.scanerrorthreshold=10;
  param.checkimages.scanerrorpos=pipeline_checkforscanerrors(param);
  save(parameterfile,'param');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write out a text file to report the scan faults
if flag.writescanfaultfile==1
  disp('---------------------------------------------------------------------');
  disp('Writing out a text file that reports scan faults...');
  disp('---------------------------------------------------------------------');
  param.writescanfaultfile.filename='hip29_scanfaults.txt';
  pipeline_writescanfaultfile(param);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image downscaling
if flag.downscale==1
  disp('---------------------------------------------------------------------');
  disp('Downscaling images ...');
  disp('---------------------------------------------------------------------');
  param.downscale.scale=16; %scaled images are downscaled to 1/8
  param.downscale.invert=0; %whether scaled images are to be inverted or not
  param.downscale.targetfileformat='png';
  [param.rawsize,param.scaledsize]=pipeline_resize(param); %writes out downscaled images to the scaleddir directory
  save(parameterfile,'param');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute absolute orientation of individual tiles (whole-tile) from cutting marks
if flag.absrot==1
  disp('---------------------------------------------------------------------');
  disp('Computing absolute orientation from cutting marks ...');
  disp('---------------------------------------------------------------------');
  param.absrot.lfreq=30/255; %0.21; %lowest frequency to consider, as percentage of image frequencies [0..1]
  param.absrot.hfreq=150/255; %0.35; %highest frequency to consider, as percentage of image frequencies [0..1]
  param.absrot.stopeach=0; %If this is 1, then the function will wait for "return" after each slice
  param=pipeline_getabsrot_fromcuts_fast_fft(param);
  param.absrotang=param.absrot.uncorrang;
  save(parameterfile,'param');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correct absolute orientation estimates BETWEEN-SLICE
if flag.correctabsrotbetween==1
  disp('---------------------------------------------------------------------');
  disp('Correcting absolute orientation estimates between-slice ...');
  disp('---------------------------------------------------------------------');
  param.correctabsrotbetween.lfreq=30/255; %0.21; %lowest frequency to consider, as percentage of image frequencies [0..1]
  param.correctabsrotbetween.hfreq=150/255; %0.35; %highest frequency to consider, as percentage of image frequencies [0..1]
  param.correctabsrotbetween.outlierthreshold=4; %3 means 3 times STD
  param=pipeline_correctabsrotbetween(param);
  param.absrotang=param.correctabsrotbetween.ang;
  save(parameterfile,'param');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correct absolute orientation estimates WITHIN-SLICE
if flag.correctabsrotwithin==1
  disp('---------------------------------------------------------------------');
  disp('Correcting absolute orientation estimates within-slice ...');
  disp('---------------------------------------------------------------------');
  param.correctabsrotwithin.in_absrotang=param.correctabsrotbetween.ang;
  param.correctabsrotwithin.in_absrotpwr=param.correctabsrotbetween.pwr;
  param.correctabsrotwithin.outlierthreshold=1; %5 means within-slice STD is larger than 5 degrees
  param.correctabsrotwithin.recomputealltiles=0; %1: outliers are detected by high STD, all tiles in that slice will be recomputed. 0: tiles which are too far off the mean for that slice will be recomputed 
  param.correctabsrotwithin.stopeach=0; %If this is 1, then the function will wait for "return" after each tile
  param=pipeline_correctabsrotwithin(param);
  param.absrotang=param.correctabsrotwithin.absrotang;
  save(parameterfile,'param');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute same-tile between-slice translations based on absrot
if flag.abstransbetweenslices==1     %Compute relative rigid translations based on absrot
  disp('---------------------------------------------------------------------');
  disp('Computing same-tile between-slice translations based on absrot ...');
  disp('---------------------------------------------------------------------');
  %param.abstransbetweenslices.in_absrotang=param.correctabsrotbetween.ang;
  %param.abstransbetweenslices.in_absrotpwr=param.correctabsrotbetween.pwr;
  param.abstransbetweenslices.in_absrotang=param.absrotang;
  param.abstransbetweenslices.dofiltering=1;
  param.abstransbetweenslices.lowestpixpercyc=0; %highest frequency
  param.abstransbetweenslices.highestpixpercyc=100; %lowest frequency
  param=pipeline_gettransbetweenslicesfromabsrot(param);
  save(parameterfile,'param');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute between-tile same-slice translations based on absrot
if flag.abstransneighbors==1     %Compute relative rigid translations based on absrot
  disp('---------------------------------------------------------------------');
  disp('Computing same-slice (neighbor tile) translations based on absrot ...');
  disp('---------------------------------------------------------------------');
  param.abstransneighbors.lowestpixpercyc=0; %highest frequency when filtering
  param.abstransneighbors.highestpixpercyc=5; %lowest frequency when filtering
  param.abstransneighbors.docomputation=0;
  param.abstransneighbors.consistencythreshold=5;
  param=pipeline_getneighbortransabsrot(param);
  save(parameterfile,'param');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute relative orientations from slice to slice
if flag.relrot==1
  disp('---------------------------------------------------------------------');
  disp('Computing relative slice orientations and translations ...');
  disp('---------------------------------------------------------------------');
  param.relrot.scale=8; %16; %downscaling factor during first rough alignment
  param.relrot.nrofangles=180; %nr of angles to be tested through 360deg
  param.relrot.dofiltering=1; %Whether images should be filtered (before further downscaling)
  param.relrot.lowestpixpercyc=0; %highest frequency when filtering
  param.relrot.highestpixpercyc=50; %100; %lowest frequency when filtering
  param.relrot.stopeach=0;
  param=pipeline_getrelrot(param);
  save(parameterfile,'param');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check relative-rotation-consistency, if there is more than one tile
if flag.relrotcheckconsistency==1
  disp('---------------------------------------------------------------------');
  disp('Checking consistency of relative slice orientations and translations ...');
  disp('---------------------------------------------------------------------');
  param.checkrelrotconsistency.reestimatestdthres=5; %re-estimate for std > 5 degrees
  param.checkrelrotconsistency.scale=2; %16; %downscaling factor during first rough alignment
  param.checkrelrotconsistency.nrofangles=180; %nr of angles to be tested through 360deg
  param.checkrelrotconsistency.dofiltering=1; %Whether images should be filtered (before further downscaling)
  param.checkrelrotconsistency.lowestpixpercyc=0; %highest frequency when filtering
  param.checkrelrotconsistency.highestpixpercyc=20; %100; %lowest frequency when filtering
  param.checkrelrotconsistency.stopeach=0;
  param=pipeline_checkrelrotconsistency(param);
  save(parameterfile,'param');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refine relative orientations from slice to slice
if flag.refinerelrot==1
  disp('---------------------------------------------------------------------');
  disp('Refining relative slice orientations and translations ...');
  disp('---------------------------------------------------------------------');
  param.refinerelrot.scale=4; %downscaling factor during first rough alignment
  param.refinerelrot.range=10; %use local range of 10 degrees (+-5)
  param.refinerelrot.nrofangles=100; %nr of angles to be tested through the range
  param.refinerelrot.stopeach=0;
  param.refinerelrot.corotarr=param.checkrelrotconsistency.rotarr;  %use corrected estimates for coarse rot and trans
  %param.refinerelrot.cotransxarr=param.checkrelrotconsistency.transxarr;
  %param.refinerelrot.cotransyarr=param.checkrelrotconsistency.transyarr;
  param=pipeline_refinerelrot(param);
  save(parameterfile,'param');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Review Rigid alignment
if flag.reviewrigid==1
  disp('---------------------------------------------------------------------');
  disp('Reviewing rigid slice orientations and translations ...');
  disp('---------------------------------------------------------------------');
  %This uses param.absrotang for absolute rotations
  pipeline_reviewrigid(param);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute relative and absolute rigid transformation matrices from absolute or relative estimates
if flag.computerigidmatrices==1
  disp('---------------------------------------------------------------------');
  disp('Compute relative and absolute rigid transformation matrices ...');
  disp('---------------------------------------------------------------------');
  param.computerigidmatrices.userelative=0; %if this is 1, relative rotation/translation estimates are used, if 0, absolute
  %param.computerigidmatrices.checkwithinslicetilealignment=0;  %This only works if there are at least 2 rows and 2 columns
  %param.computerigidmatrices.withinslicethreshold=10; %if the mean STD for predicted tile corner positions is above this (in downscaled pixels), position is corrected
  %param.computerigidmatrices.checkbetweenslicetilealignment=1;
  param.computerigidmatrices.breakthreshold=100; %if the distance between predictions is larger than this (in downscaled pixels), assume a break
  param=pipeline_computerigidmatrices(param);  %This uses param.absrotang for absolute rotations
  save(parameterfile,'param');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Render Rigid alignment, based on the matrices computed by computerigidmatrices
if flag.renderrigid==1
  disp('---------------------------------------------------------------------');
  disp('Rendering a stack of rigidly aligned images ...');
  disp('---------------------------------------------------------------------');
  param.renderrigid.targetdir='../../rigid/';
  %param.renderrigid.targetdir='../../rigidblended/';
  param.renderrigid.basename='W2W7_rigid16_s%03d.png';
    %rendermode 1: single tile per image, 
    %rendermode 2: all tiles of one slice in one image, colored, 
    %rendermode 3: all tiles of one slice in one image, blended
  param.renderrigid.rendermode=1; 
  param.renderrigid.startslice=1;
  param.renderrigid.endslice=param.nrofslices;
    
  pipeline_renderrigid(param);
  save(parameterfile,'param');
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Predict where corresponding points should be on neighboring slices, based on rigid alignment
if flag.predictcorr==1
  disp('---------------------------------------------------------------------');
  disp('Predict corresponding point positions based on rigid alignment ...');
  disp('---------------------------------------------------------------------');
  param.predictcorr.gridwidth=8;
  param.predictcorr.gridheight=8;
  param.predictcorr.patchwidth=param.rawsize(2)/param.predictcorr.gridwidth; %this will cover the whole src image
  param.predictcorr.patchheight=param.rawsize(1)/param.predictcorr.gridheight; %CAUTION - in raw-scale!
  param=pipeline_predictcorr(param);
  save(parameterfile,'param');
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute actual peak of cross-correlation between local image regions
% Uses the location estimates computed by pipeline_predictcorr from the rigid alignment
% For regions close to the edge, the algorithm attempts to move the region inside the image
if flag.computelocalcorr==1
  disp('---------------------------------------------------------------------');
  disp('Compute cross-correlation for local regions based on rigid alignment ...');
  disp('---------------------------------------------------------------------');
  %Gridwidth,gridheight,patchwidth and patchheight are taken from predictcorr
  param.computelocalcorr.minoverlap=0.1; %between 0 and 1 (relative to patch size; patches on image borders are used if the overlap is larger than minoverlap)
  param.computelocalcorr.showimages=1;
  param.computelocalcorr.stopeach=0;
  param.computelocalcorr.dofiltering=1;
  param.computelocalcorr.minppcyc=1;
  param.computelocalcorr.maxppcyc=200;
  param.computelocalcorr.downscalefactor=8; %2 means divided by 2
  param=pipeline_computelocalcorr(param);
  save(parameterfile,'param');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit AFFINE transformations to the local correspondences iteratively
if flag.fititerativeaffine==1
  disp('---------------------------------------------------------------------');
  disp(' Fitting AFFINE transformations to the local correspondences iteratively...');
  disp('---------------------------------------------------------------------');
  param.fititerativeaffine.vargt075threshold=100; %Corresponding points are only considered if their vargt075 value is below this threshold
  param.fititerativeaffine.showprealignment=0;
  param.fititerativeaffine.nrofiterations=100;
  param=pipeline_fititerativeaffine(param);
  save(parameterfile,'param');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Render stack of AFFINE aligned images
if flag.renderaffine==1
  disp('---------------------------------------------------------------------');
  disp(' Rendering out a stack of affinely aligned images...');
  disp('---------------------------------------------------------------------');
  %param.renderaffine.targetdir='../../rigidcolor/';
  param.renderaffine.targetdir='../../affine_big/';
  param.renderaffine.basename='W2W7_affine_s%03d.png';
  param.renderaffine.usedownscaled=0;
    %rendermode 1: single tile per image, 
    %rendermode 2: all tiles of one slice in one image, colored, 
    %rendermode 3: all tiles of one slice in one image, blended
  param.renderaffine.rendermode=1; 
  param.renderaffine.startslice=183;
  param.renderaffine.endslice=param.nrofslices;
  param.renderaffine.usetiledrendering=1;
  param.renderaffine.tilesize=1024;
  
  if param.renderaffine.usetiledrendering==1
    pipeline_renderaffinetiled(param);
  else
    pipeline_renderaffine(param);
  end;
  save(parameterfile,'param');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Render cropped region of stack of AFFINE aligned images
if flag.renderaffineregion==1
  disp('---------------------------------------------------------------------');
  disp(' Rendering out a local region of a stack of affinely aligned images...');
  disp('---------------------------------------------------------------------');
  %param.renderaffine.targetdir='../../rigidcolor/';
  param.renderaffineregion.targetdir='../../affineregion/';
  param.renderaffineregion.basename='W2W7_affinecropped_s%03d.png';
  param.renderaffineregion.usedownscaled=0;
    %rendermode 1: single tile per image, 
    %rendermode 2: all tiles of one slice in one image, colored, 
    %rendermode 3: all tiles of one slice in one image, blended
  %param.renderaffine.rendermode=1; 
  param.renderaffineregion.startslice=1;
  param.renderaffineregion.endslice=param.nrofslices;
  param.renderaffineregion.x1=5000;
  param.renderaffineregion.y1=5000;
  param.renderaffineregion.width=512;
  param.renderaffineregion.height=512;
  pipeline_renderaffineregion(param);
  save(parameterfile,'param');
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute actual peak of cross-correlation between local image regions
% Uses the location estimates computed by pipeline_fititerativeaffine from the affine alignment
% For regions close to the edge, the algorithm attempts to move the region inside the image
if flag.computefinelocalcorr==1
  disp('---------------------------------------------------------------------');
  disp('Compute cross-correlation for fine local regions based on affine alignment ...');
  disp('---------------------------------------------------------------------');
  %Gridwidth,gridheight,patchwidth and patchheight are taken from predictcorr
  param.computefinelocalcorr.gridwidth=64;
  param.computefinelocalcorr.gridheight=64;
  param.computefinelocalcorr.patchwidth=param.rawsize(2)/param.computefinelocalcorr.gridwidth; %this will cover the whole src image
  param.computefinelocalcorr.patchheight=param.rawsize(1)/param.computefinelocalcorr.gridheight; %CAUTION - in raw-scale!
  param.computefinelocalcorr.minoverlap=0.1; %between 0 and 1 (relative to patch size; patches on image borders are used if the overlap is larger than minoverlap)
  param.computefinelocalcorr.showimages=0;
  param.computefinelocalcorr.stopeach=0;
  param.computefinelocalcorr.dofiltering=0;
  param.computefinelocalcorr.minppcyc=1;
  param.computefinelocalcorr.maxppcyc=200;
  param.computefinelocalcorr.downscalefactor=4; %2 means divided by 2
  [param,finecorrparam]=pipeline_computefinelocalcorr(param);
  save(parameterfile,'param');
  save(finecorrparamfile,'finecorrparam','-V7.3');
end;

