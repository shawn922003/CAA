function [num, den] = pade_from_moments_sym(mrow, L, M)
% mrow: [m0 m1 ... m_{L+M}]（row 向量, double 或 sym 都可）
% 回傳：
%   num(1,L+1), den(1,M+1)  使得  Pade(s) = (num0+num1*s+...+numL*s^L) / (1+den1*s+...+denM*s^M)

    K = L + M;
    mrow = mrow(:).';                          % 保證 row
    if numel(mrow) < K+1
        error('Need at least L+M+1 moments');
    end
    m0 = sym(0);                               % 確保是 sym 0

    % ---- 解分母 a1..aM：
    % 對 k = L+1..L+M，滿足 sum_{i=1..M} a_i * m_{k-i} = - m_k
    % 其中當 k-i < 0 時，m_{k-i} 視為 0
    A = sym(zeros(M, M)); 
    b = sym(zeros(M, 1));
    for r = 1:M
        k = L + r;                              % 這裡的 k 是「階數」
        for c = 1:M
            idx = k - c + 1;                    % 對應 m_{k-c} 的 1-based 索引
            if idx >= 1
                A(r,c) = mrow(idx);
            else
                A(r,c) = m0;                    % 越界→0
            end
        end
        b(r) = -mrow(k + 1);                    % -m_k
    end
    a = A \ b;                                  % a = [a1..aM].'
    den = [sym(1), a.'];                        % Q(s) = 1 + a1 s + ... + aM s^M

    % ---- 分子 b_k（k=0..L）：
    % b_k = sum_{i=0}^k m_i * a_{k-i}，其中 a_0=1、a_j=0 (j>M)
    num = sym(zeros(1, L+1));
    for k = 0:L
        acc = m0;
        for i = 0:k
            j = k - i;                          % 對應 a_j
            if j == 0
                aj = sym(1);
            elseif j <= M
                aj = a(j);
            else
                aj = m0;                        % 超過 M → 0
            end
            acc = acc + mrow(i + 1) * aj;       % m_i * a_{k-i}
        end
        num(k+1) = acc;
    end
end
