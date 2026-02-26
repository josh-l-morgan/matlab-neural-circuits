function[num] = getFirstNumber(nam)

 %%Get cell id numbers
   isStart = 0;
   ns = char;
   for s = 1:length(nam)
      n = str2num(nam(s));
      if isStart & isempty(n)
          break
      end
      if ~isempty(n) 
      
          isStart = 1;
          ns = [ns nam(s)];
      end
   end
   
   num = str2num(ns);
   if isempty(num)
       num = 0;
   end
