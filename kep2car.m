function [r,v] =kep2car(kepl_coord,mu)

%it gives r and v in keplerian frame as vectors
% INPUT: kep_coord = [a,e,i,omega,w,f]: i,omega,w,f must be in rad
%OUTPUT: r and v are column vectors

a = kepl_coord(1);
e = kepl_coord(2);
i = kepl_coord(3);
omega = kepl_coord(4);
w = kepl_coord(5);
f = kepl_coord(6);

h = sqrt(a*(1-e^2)*mu);
r_peri = h^2/(mu*(1+e*cos(f)))*[cos(f); sin(f); 0];
v_peri = mu/h*[-sin(f); e+cos(f); 0];

R3_w = [cos(w), sin(w), 0; -sin(w), cos(w),0; 0,0,1];
R1_i = [1,0,0; 0, cos(i), sin(i); 0, -sin(i), cos(i)];
R3_omega = [cos(omega), sin(omega),0; -sin(omega), cos(omega), 0; 0,0,1];

r = (R3_w*R1_i*R3_omega)'*r_peri;
v = (R3_w*R1_i*R3_omega)'*v_peri;


