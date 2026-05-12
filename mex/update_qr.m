function [Rt, Qhat2] = update_qr(Hcheck, Qcheck2)
%UPDATE_QR Updates the QR factorisation of the condensed reduced order MPC
%problem.
%   Hcheck:     Upper Hessenberg matrix
%               (N*m+1 by N*m) matrix
%   Qcheck:     Part of the QR decomposition of E'