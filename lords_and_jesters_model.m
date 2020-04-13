% Set constant parameters
a = .8;
kappa = a*((1 - a)/a)^(1 - a);
lkt = .05;
lft = .95;
b = .2;
B = b^b*(1 - b)^(1 - b);
%% Farmers as jesters
% Loop through different values of subsistence level
m=(0:50)/100;
X_f = zeros(length(m),2);
for i=1:length(m)
    sub = m(i);
    X_f(i,:) = fminsearch(@(w) getWage(w,a,b,B,sub,lkt,lft,kappa),[.5,1]);
end
%% Kings as jesters
% Loop through different values of subsistence level
X_k = zeros(length(m),2);
for i=1:length(m)
    sub = m(i);
    W = fminsearch(@(w) getWage(w,a,1-b,B,sub,lft,lkt,kappa),[.5,1]);
    X_k(i,1) = W(2);
    X_k(i,2) = W(1);    
    % Note that the wages returned are now in the opposite order 
end
% Calculate total labor by type
lt_f = X_f.*X_f(:,1).^(a-1).*kappa;
lt_k = X_k.*X_k(:,2).^(a-1).*kappa;
% Calculate consumption of corn
ct_f = a.*lt_f.*X_f+(1-a).*m';
ct_k = a.*lt_k.*X_k+(1-a).*m';
%Calculate consumption of jesters
jt_f = ((1-a).*lt_f.*X_f-(1-a).*m')./X_f(:,1);
jt_k = ((1-a).*lt_k.*X_k-(1-a).*m')./X_k(:,2);
%Calculate utility
ut_f = (ct_f-m').^a.*jt_f.^(1-a)-1./2.*lt_f.^2;
ut_k = (ct_k-m').^a.*jt_f.^(1-a)-1./2.*lt_k.^2;

figure

plot(m,X_f(:,1),m,X_f(:,2),m,X_k(:,1),m,X_k(:,2))
axis([0 .5 0 1.5])
legend('Farm Jesters-Farm Wage','Farm-Jesters-King Wage',...
'King Jesters-Farm Wage','King-Jesters-King Wage')

figure

plot(m,lt_f(:,1),m,lt_f(:,2),m,lt_k(:,1),m,lt_k(:,2))
axis([0 .5 0 1.5])
legend('Farm Jesters-Farm Effort','Farm-Jesters-King Effort',...
'King Jesters-Farm Effort','King-Jesters-King Effort')

figure 
plot(m,ct_f(:,1)./(X_f(:,1).*lt_f(:,1)),...
    m,ct_f(:,2)./(X_f(:,2).*lt_f(:,2)),...
    m,ct_k(:,1)./(X_k(:,1).*lt_k(:,1)),...
    m,ct_k(:,2)./(X_k(:,2).*lt_k(:,2)))
axis([0 .5 0 1.5])
legend('Farm Jesters-Farm Corn','Farm-Jesters-King Corn',...
'King Jesters-Farm Corn','King-Jesters-King Corn')
% Utility
figure 
plot(m,ut_f(:,1),...
    m,ut_f(:,2),...
    m,ut_k(:,1),...
    m,ut_k(:,2))
legend('Farm Jesters-Farm Utility','Farm-Jesters-King Utility',...
'King Jesters-Farm Utility','King-Jesters-King Utility')
%% Incorporate taxes by income with farmers as jesters
% Loop through different values of subsistence level and taxes
t = (0:90)/100;
m=(0:10)/100;
Xf_f = zeros(length(m),length(t));
Xk_f = zeros(length(m),length(t));
for i=1:length(m)
    for j = 1: length(t)
        sub = m(i);
        tax = t(j);
        W = fminsearch(@(w) getWageWithTaxesFarmJesters(w,a,b,B,sub,lkt,lft,kappa,tax),[.5,1]);
        Xf_f(i,j) = W(1);
        Xk_f(i,j) = W(2);
    end
end

% Calculate price of jesters
p = Xf_f;
% Calculate total labor by type
lf_f = (Xf_f.*p.^(a-1).*kappa).^.5;
lk_f = (Xk_f.*(1-t).*p.^(a-1).*kappa).^.5;
% Calculate consumption of corn
cf_f = a.*lf_f.*Xf_f+a.*Xk_f.*lk_f.*t*lkt/lft+(1-a).*m';
ck_f = a.*lk_f.*Xk_f.*(1-t)+(1-a).*m';
%Calculate consumption of jesters
jf_f = ((1-a).*lf_f.*Xf_f+(1-a).*Xk_f.*lk_f.*t*lkt/lft-(1-a).*m')./p;
jk_f = ((1-a).*lk_f.*Xk_f.*(1-t)-(1-a).*m')./p;
%Calculate utility
uf_f = (cf_f-m').^a.*jf_f.^(1-a)-1./2.*lf_f.^2;
uk_f = (ck_f-m').^a.*jk_f.^(1-a)-1./2.*lk_f.^2;
% Get Total Utility
U = uf_f.*lft+uk_f.*lkt;
% Plot utility
figure
plot(t,U(1,:))
figure
plot(t,U(length(m),:))

%% Incorporate taxes by income with kings as jesters
t = (0:90)/100;
Xf_k = zeros(length(m),length(t));
Xk_k = zeros(length(m),length(t));
for i=1:length(m)
    for j = 1: length(t)
        sub = m(i);
        tax = t(j);
        W = fminsearch(@(w) getWageWithTaxesKingJesters(w,a,b,B,sub,lkt,lft,kappa,tax),[.5,1]);
        Xf_k(i,j) = W(1);
        Xk_k(i,j) = W(2);
    end
end

% Calculate price of jesters
p = Xk_k;
% Calculate total labor by type
lf_k = (Xf_k.*p.^(a-1).*kappa).^.5;
lk_k = (Xk_k.*(1-t).*p.^(a-1).*kappa).^.5;
% Calculate consumption of corn
cf_k = a.*lf_k.*Xf_k+a.*Xk_k.*lk_k.*t*lkt/lft+(1-a).*m';
ck_k = a.*lk_k.*Xk_k.*(1-t)+(1-a).*m';
%Calculate consumption of jesters
jf_k = ((1-a).*lf_k.*Xf_k+(1-a).*Xk_k.*lk_k.*t*lkt/lft-(1-a).*m')./p;
jk_k = ((1-a).*lk_k.*Xk_k.*(1-t)-(1-a).*m')./p;
%Calculate utility
uf_k = (cf_k-m').^a.*jf_k.^(1-a)-1./2.*lf_k.^2;
uk_k = (ck_k-m').^a.*jk_k.^(1-a)-1./2.*lk_k.^2;
% Get Total Utility
U = uf_k.*lft+uk_k.*lkt;
% Plot utility
figure
plot(t,U(1,:))
figure
plot(t,U(length(m),:))
%% Kings as Jesters Taxing Jester 
wf_jester_tax = zeros(length(m),length(t));
wk_jester_tax = zeros(length(m),length(t));
J = zeros(length(m),length(t));
for i=1:length(m)
    for j = 1: length(t)
        sub = m(i);
        tax = t(j);
        W = fminsearch(@(w) getWageWithJesterTaxKingJesters(w,a,b,B,sub,lkt,lft,kappa,tax),[.5,1,.5]);
        wf_jester_tax(i,j) = W(1);
        wk_jester_tax(i,j) = W(2);
        J(i,j) = W(3);
    end
end

% Calculate price of jesters
p = (1+t).*wk_jester_tax;
% Calculate amount of labor 
lf = (kappa.*wf_jester_tax.*p.^(a-1)).^(1/2);
lk = (kappa.*wk_jester_tax.*p.^(a-1)).^(1/2);
%Calculate consumption of jesters
jf = ((1-a).*(lf.*wf_jester_tax+t./(1+t).*p.*J/lft)-(1-a).*m')./p;
jk = ((1-a).*lk.*wk_jester_tax-(1-a).*m')./p;    
% Calculate consumption of corn
cf = a.*(lf.*wf_jester_tax+t./(1+t).*p.*J/lft)+(1-a).*m';
ck = a.*lk.*wk_jester_tax+(1-a).*m';    
%Calculate utility
uf_k = (cf-m').^a.*jf.^(1-a)-1./2.*lf.^2;
uk_k = (ck-m').^a.*jk.^(1-a)-1./2.*lk.^2;
% Get Total Utility
U2 = uf_k.*lft+uk_k.*lkt;
% Plot utility
figure
plot(t,U2(1,:))
figure
plot(t,U2(length(m),:))