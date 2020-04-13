function error = getWageWithTaxesKingJesters(w,a,b,B,m,lkt,lft,kappa,t)
    wf0 = w(1);
    wk0 = w(2);
    % Calculate jester price
    p = wk0;
    % Calculate amount of labor 
    lf = (kappa*wf0*p^(a-1))^(1/2);
    lk = (kappa*wk0*(1-t)*p^(a-1))^(1/2);
    % Calculate consumption of corn
    cf = a*(lf*wf0+wk0*lk*t*lkt/lft)+(1-a).*m;
    ck = a*(1-t)*lk*wk0+(1-a)*m;       
    
    error = [B*wk0^(-b)-wf0^(1-b);...
        lft*lf*wf0/(1-b)-lkt*ck-lft*cf];
    error = error'*error;
end