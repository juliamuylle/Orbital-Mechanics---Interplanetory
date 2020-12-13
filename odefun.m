function dy = odefun (nu,y,t)

r = y(1:3);
dy = [y(4:6); -nu/(norm(r)^3)*y(1:3)];

end