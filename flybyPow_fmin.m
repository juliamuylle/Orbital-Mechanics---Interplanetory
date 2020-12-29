function [deltav_perig,rp] = flybyPow(vinfM,vinfP,mu_E,hatm,Re)
%[DV_ga,DV_fb,rp,steps,delta,arcs]

delta_true = acos(dot(vinfM,vinfP)/(norm(vinfM)*norm(vinfP)));

deta_deg = rad2deg(delta_true);

eM = @(rp) 1+((rp*(norm(vinfM))^2)/mu_E);
deltaM = @(rp) 2*asin(1./eM(rp));

eP = @(rp) 1+((rp*(norm(vinfP))^2)/mu_E);
deltaP = @(rp) 2*asin(1./eP(rp));

delta = @(rp) deltaM(rp)/2+deltaP(rp)/2;

FUN = @(rp) delta(rp)-delta_true;

% x = linspace(0,100000,100000);
% figure(30)
% plot(x,FUN(x))
% hold on

options = optimset('TolFun',1e-14,'Display','off');
rp = fsolve(FUN,Re+hatm,options);

%                 if rp <= Re+hatm
%                     rp = NaN;
%                 end

eM_val = eM(rp);
eP_val = eP(rp);

hM = sqrt(mu_E*rp*(1+eM_val));
vperig_M = mu_E/hM*(1+eM_val);

hP = sqrt(mu_E*rp*(1+eP_val));
vperig_P = mu_E/hP*(1+eP_val);

deltav_perig = abs(vperig_P-vperig_M);


