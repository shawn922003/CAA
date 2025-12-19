function [Q, H] = arnoldi(A, q1, m)
%ARNOLDI_SIMPLE  Minimal Arnoldi (as in the provided pseudocode).
%
% Inputs:
%   A  : n-by-n matrix
%   q1 : n-by-1 starting vector (will be normalized)
%   m  : number of Arnoldi vectors to generate (Q will be n-by-m)
%
% Outputs:
%   Q  : n-by-m orthonormal basis (q1..qm)
%   H  : m-by-(m-1) Hessenberg entries h_{j,k-1}

    q1 = q1(:);
    q1 = q1 / norm(q1);

    n = length(q1);
    Q = zeros(n, m);
    H = zeros(m, m-1);

    Q(:,1) = q1;

    for k = 2:m
        qk = A * Q(:,k-1);

        for j = 1:(k-1)
            H(j, k-1) = Q(:,j)' * qk;   % h_{j,k-1}
            qk = qk - H(j, k-1) * Q(:,j);
        end

        H(k, k-1) = norm(qk);           % h_{k,k-1}

        % 若發生 breakdown（norm 太小），這裡簡單報錯或停止
        if H(k, k-1) == 0
            error('Arnoldi breakdown at k = %d (hk,k-1 = 0).', k);
        end

        Q(:,k) = qk / H(k, k-1);
    end
end
