function error = getWage(w,a,b,B,m,lkt,lft,kappa)
    wf0 = w(1);
    wk0 = w(2);
    error = [B*wk0^(-b)-wf0^(1-b);
        lkt*kappa*(1/b-a)*wk0^2*wf0^(a-1)-...
        lft*kappa*a*wf0^(a+1)-...
        (lkt+lft)*(1-a)*m];
    error = error'*error;
end