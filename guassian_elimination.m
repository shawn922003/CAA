function [x, status] = guassian_elimination(A, b, tol)
% Solve A x = b by Gaussian elimination with partial pivoting
% status: "ok" | "singular"
    if nargin < 3, tol = 1e-12; end

    [m, n] = size(A);
    b = b(:);  % 保證是 column vector
    if m ~= n || n ~= numel(b)
        error('Dimension mismatch: A must be n×n and b must be n×1');
    end

    % --- Forward elimination with partial pivoting ---
    for i = 1:n
        % 1) pivot: 選 |A(k,i)| 最大者與第 i 列互換
        [~, p_rel] = max(abs(A(i:m, i)));
        p = p_rel + i - 1;
        if abs(A(p, i)) < tol
            x = NaN(n,1);
            status = "singular";
            return
        end
        if p ~= i
            A([i p], :) = A([p i], :);
            b([i p])   = b([p i]);
        end

        % 2) 消去（向量化）
        rows = (i+1):n;
        if ~isempty(rows)
            factors = A(rows, i) / A(i, i);          % (n-i)×1
            A(rows, i:n) = A(rows, i:n) - factors * A(i, i:n);
            A(rows, i)   = 0;                        % 清掉數值殘留
            b(rows)      = b(rows)      - factors * b(i);
        end
    end

    % --- Back substitution ---
    x = zeros(n, 1);
    for i = n:-1:1
        s = A(i, i+1:n) * x(i+1:n);
        x(i) = (b(i) - s) / A(i, i);
    end
    status = "ok";
end
