function [c,s] = givens(a,b) %# codegen
%GIVENS Returns the Given rotation for the two scalars a and b, such that
%R^\top [a;b] = [c; 0], with R = [c s; -s c]. It follows [Golub, Alg.
%5.1.3] 
%   a,b:    Real numbers
%   c,s:    Real numbers

if b==0
    c=1; s=0;
else
    if abs(b) > abs(a)
        t = -a/b;
        s = 1 / sqrt(1+t^2);
        c = s*t;
    else
        t = -b/a;
        c = 1 / sqrt(1+t^2);
        s = c*t;
    end
end
end