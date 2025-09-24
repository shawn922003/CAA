function C = taylor_coeffs(f, x, a, K)
% 取 f 在 x=a 的泰勒係數：回傳 1×(K+1) 向量 C，C(n+1)=c_n
% f: sym（可為純量或向量）；x: sym 變數；a: 展開點；K: 階數
    f = sym(f);
    sz = size(f);
    C = sym(zeros(numel(f), K+1));
    for i = 1:numel(f)
        fi = f(i);
        for n = 0:K
            C(i, n+1) = subs(diff(fi, x, n), x, a) / factorial(n);
        end
    end
    C = reshape(C, [sz(1), sz(2)*(K+1)]);  % 方便看也可不 reshape
end