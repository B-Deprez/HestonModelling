function [C, Ceq] = mycon(k, n, t, v, r)
% This function is used to express the inequality constraint C <= 0 for the
% optimisation
    C = t^2-2*k*n + 0*v + 0*r;
    Ceq = [];
    