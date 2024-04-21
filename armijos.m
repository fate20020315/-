

%%%%%%%%%%%%%%%%%%%%%%%%%%%armjioÏßËÑË÷%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function skk= armijos(fun, xk,k, rho, d)  
c=0.3;
% sk=1/sqrt(1-0.01*k);
mk = 0; 
x = xk - d;
if fun(x) <=  (1-c*rho^(mk))*fun(xk )+1/((k)^(4/3))
        sk=1;
    else
       sk=1/sqrt(k);
     % sk=1/k; 
      %sk=0.3;
       
end
skk=sk;
end



