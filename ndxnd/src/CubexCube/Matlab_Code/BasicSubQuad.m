function[s,ws]= BasicSubQuad(j,k,d,t,wt,N)   
     %Step 5: Transform to pyramids
      [s,ws]=Step5(j,k,d,t,wt);
     %Step 4: Reflections
      [s,ws]=Step4(N,k,s,ws);
end
