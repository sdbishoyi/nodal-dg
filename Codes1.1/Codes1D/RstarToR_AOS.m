function [Raos,RaosrH] = RstarToR_AOS(x,rH,e)

%Uses newton's method to turn tortoise coordinates to standard schwarschild
%coordinates.
%rstar = read in rstar coordinates, raos = AOS Black Hole r coordinate..but
%computed iteratively using newtons method.
%raosrH + rH ~ raos... so raosrH near 0 at horizon instead of rH,
%we actually update raosrH and add rH to this number.
%rsn = rstar (newton).
%parameter 'e' and 'L' stands for epsilon, and capital Lambda respectively.
% Most things have been commented out to check functionality.

      Raos = zeros(size(x));
      RaosrH = zeros(size(x));
      LENGTHx = size(x,1)*size(x,2);
      
      for k = 1:LENGTHx  %finds rsch at each point
          rstar = x(k);
          if  rstar < -1.0e3   %less than -1000, very close to horizon (when m=1)
              error = 'Tortoise coordinate too small'; %#ok<NASGU> 
              pause
          end

        if  rstar > rH           %not near horizon
            raos = rstar;
            raosrH = raos - rH;        
        else                     %near horizon                           
             raosrH = (rH/(1+e))*exp(rstar/((1+e)^2*rH) - 1);              
%            raosrH = rH*((1 + exp(rstar/((1+e)^2*rH) - 1))^(1/(1+e))-1);
             raos   = raosrH + rH;                           
        end
        for i=1:100              %Newtons method
        
         hgfunc =  real(hypergeom([1 2/(1+e)],1 + 2/(1+e),(raos/rH).^(1+e)));                    
         rsn  = -0.5*(1+e)*(raos.^2/rH).*hgfunc;  
         
         if abs(rsn-rstar) < 1.0e-25  %approximation is good enough
              break             
         end                              
         c1 = 2*raos*real(hypergeom([1 2/(e + 1)], 2/(e + 1) + 1, (raos/rH).^(e + 1)));
         d1 = real(hypergeom([2 2/(e + 1) + 1], 2/(e + 1) + 2, (raos/rH).^(e + 1)));
%        d2 = rH*(2/(e + 1) + 1);
         c2 = 2*raos^2*(raos/rH)^e*(1+e)/(3+e)*d1/rH; %/d2;
         drsdr = -((e + 1)*(c1 + c2))/(2*rH); 
         
               rtemprH = (raos-rH) + (rstar-rsn)/drsdr;      
               if rtemprH < 0
                  rtemprH = 0.5*(raos - rH);
               end
             raosrH  = rtemprH;
             raos     = raosrH + rH;
        end
       Raos(k) = raos;
       RaosrH(k) = raosrH;
       
      end