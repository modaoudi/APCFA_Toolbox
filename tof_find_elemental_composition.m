function [elemComp, Isel] = tof_find_elemental_composition(m,tol, mss)
%
%
%

Isel = mss(:,1) <= (m+tol) & mss(:,1) >= (m-tol);

elemComp = mss(Isel,1:end);