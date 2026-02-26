
sourceFileName = 'F:\LindaXu\rayAlign01\Segmentation1-LX.vss'
newName = 'Segmentation1-LX'
TPN = 'L:\joshm\Segmentation\Vast\rayAlign01\LindaXu\autoBackup\'


while 1
    timeNow = clock;
   newFileName =  sprintf('%s%s_%d+%02.0f+%02.0f.vss',TPN,newName,timeNow(1),timeNow(2),timeNow(3));
   tempFileName =  sprintf('%s%s_tmp.vss',TPN,newName);
   disp(sprintf('Copying %s to %s',sourceFileName,newFileName))
   if exist(newFileName,'file')
      try copyfile(newFileName,tempFileName)
      catch err
          err
      end
   end
   try   copyfile(sourceFileName,newFileName)
   catch err
       err
   end
    
    pause(3600)
end