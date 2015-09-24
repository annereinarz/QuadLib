function N = FindSubset(numN,k)

 N=[];
 s=dec2bin(numN);
 cnt=1;
 for i=length(s):-1:1
     if(strcmp(s(i),'1'))
         N=[N cnt];
     else
         N=[N 0];
     end
     cnt=cnt+1;
 end
 if(length(s)~=k)
     N=[N zeros(1, k-length(s))];
 end

end

