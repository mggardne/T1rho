function [f,df] = exp_fun1(rp,xdata)
%EXP_FUN  Calculates a two parameter exponential function and its
% Jacobian.  For use with T1rho_map2.m.
%
% EXP_FUN(RP,XDATA) Calculates an exponential function using the two
% parameters in vector, RP, and using the independent data in vector,
% XDATA.
%
% F = EXP_FUN(RP,XDATA) Returns a vector, F, of the values of
% the exponential function for the values in XDATA.
%
% [F,DF] = EXP_FUN(RP,XDATA) Returns the Jacobian in the matrix,
% DF.
%
% NOTES:  The exponential function is:
% F(RP,XDATA) = RP(1)*EXP(-XDATA/RP(2))
%
% 10-Jul-2020 * Mack Gardner-Morse
%

%#######################################################################
%
% Calculate the Exponential Function
% F(RP,XDATA) = RP(1)*EXP(-XDATA/RP(2))
%
f = rp(1)*exp(-xdata/rp(2));
%
% Calculate the Partial Derivatives (Jacobian)
%
if nargout>1
%
  rp = rp(:);
  xdata = xdata(:);
%
  df = zeros(size(xdata,1),size(rp,1));
  df(:,1) = exp(-xdata/rp(2));
  df(:,2) = rp(1)*xdata.*exp(-xdata/rp(2))/(rp(2)*rp(2));
%
end
%
return