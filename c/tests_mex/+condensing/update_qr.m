function [Rt,Qhat2] = update_qr(Hcheck,Qcheck2) %# codegen
%UPDATE_QR Updates the QR factorisation of the condensed reduced MPC
%problem.
%   Hcheck:     Upper Hessenberg matrix
%               (N*m+1 by N*m) matrix
%   Qcheck:     Part of the QR decomposition of E'
%               (N(n+m)+s by N*m) matrix

d = size(Hcheck,2); % d = N*m

% initialise
Rt = Hcheck;
Qhat2 = [zeros(size(Qcheck2,1),1), Qcheck2; ...
    1, zeros(1,size(Qcheck2,2))];

% iterate over the entries that have to be deleted
for k=1:d
    % get Givens rotation
    [ck,sk] = condensing.givens(Rt(k,k),Rt(k+1,k));
    for j=1:size(Qhat2,1)
        % update rows k and k+1 of Rt
        if (j>=k) && (j<=d)
            tau1 = Rt(k,j); tau2 = Rt(k+1,j);
            Rt(k,j) = ck*tau1 - sk*tau2;
            Rt(k+1,j) = sk*tau1 + ck*tau2;
        end
        % update column k and k+1 of Qcheck2
        t1 = Qhat2(j,k); t2 = Qhat2(j,k+1);
        Qhat2(j,k) = ck*t1 - sk*t2;
        Qhat2(j,k+1) = sk*t1 + ck*t2;
    end
end
% delete zero row
Rt = Rt(1:end-1,:);
end