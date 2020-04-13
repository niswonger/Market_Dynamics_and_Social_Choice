function error = getWageWithJesterTaxKingJesters(w,a,b,B,m,lkt,lft,kappa,t)
    wf0 = w(1);
    wk0 = w(2);
    J = w(3);
    % Calculate jester price
    p = (1+t)*wk0;
    % Calculate amount of labor 
    lf = (kappa*wf0*p^(a-1))^(1/2);
    lk = (kappa*wk0*p^(a-1))^(1/2);
    %Calculate consumption of jesters
    jf = ((1-a)*(lf*wf0+t/(1+t)*p*J/lft)-(1-a).*m)/p;
    jk = ((1-a)*lk*wk0-(1-a)*m)/p;    
    % Calculate consumption of corn
    cf = a*(lf*wf0+t/(1+t)*p*J/lft)+(1-a).*m;
    ck = a*lk*wk0+(1-a)*m;    
    
    error1 = B*wk0^(-b)-wf0^(1-b);
    error2 = J - lkt*jk-lft*jf;
    error3 = lft*lf*wf0/(1-b)-lkt*ck-lft*cf;
    
    error = [error1;error2;error3];
    
    
    error = error'*error;
end