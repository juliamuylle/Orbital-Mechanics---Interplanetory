function [deltav_perig,rp,delta,arcs] = flybyPow(vinfM,vinfP,mu_E,rp0)
%[DV_ga,DV_fb,rp,steps,delta,arcs]

delta_true = acos(dot(vinfM,vinfP)/(norm(vinfM)*norm(vinfP)));
deta_deg = rad2deg(delta_true);

eM = @(rp) 1+(rp*(norm(vinfM))^2)/mu_E;
deltaM = @(rp) 2*asin(1/eM(rp));

eP = @(rp) 1+(rp*(norm(vinfP))^2)/mu_E;
deltaP = @(rp) 2*asin(1/eP(rp));

delta = @(rp) deltaM(rp)/2+deltaP(rp)/2;

FUN = @(rp) delta(rp)-delta_true;


rp = fzero(FUN,rp0);
%hatm = 

% if rp < Re + hatm
%     error('rp is too small');
% end

eM = eM(rp);
eP = eP(rp);

hM = sqrt(mu_E*rp*(1+eM));
vperig_M = mu_E/hM*(1+eM);

hP = sqrt(mu_E*norm(rp)*(1+eP));
vperig_P = mu_E/hP*(1+eP);

deltav_perig = abs(vperig_M-vperig_P);

delta = delta_true;

arc1 = [eM,vperig_M];
arc2 = [eP,vperig_P];
arcs = [arc1;arc2];