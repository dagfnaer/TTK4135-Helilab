function [c,ceq] = nonlconst(z)
global N;
alfa = 0.2;
betha = 20;
lamda_t = 2*pi/3;
c = [];
for i=1:6:6*N
    temp = alfa*exp(-betha*(z(i) - lamda_t)^2) - z(i+4);
    c = [c;temp];
    ceq = [];
end
end