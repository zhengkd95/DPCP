%--------------------------------------------------------------------------
% Dual Principal Component Pursuit
% see Tsakiris M C, Vidal R. Dual principal component pursuit[C]
% //Proceedings of the IEEE ICCV 2015: 10-18.
% B = DPCP(X,c,e,max_iter)
% INPUT:
% X: dimension reduced data (PCA recommended)
%    required
% c: potential rank of orthogonal space (optional)
%    default=1
% e: iteration improvement (optional)
%    default=1e-5
% maxiter: max iteration number
%    default=100
%
% OUTPUT:
% B: set of the orthogonal vectors
%--------------------------------------------------------------------------
% Copyright @ Kedi Zheng, 2018
%--------------------------------------------------------------------------
function B = DPCP(X,c,e,max_iter)

if nargin < 4
    max_iter = 1e2;
    if nargin < 3
        e = 1e-5;
        if nargin < 2
            c = 1;
        end
    end
end


B = [];
if c > 0
    k = 0; dJ = 1e3; J = [];
    [~,S,V] = svd(X');
    n0 = V(:,end);
    while (k<=max_iter && abs(dJ) > e)
        fprintf('Iter=%d, J=%g\n',k,norm(X'*n0,1))
        k = k+1;
        n = l1_min(X,n0);
        dJ = norm(X'*n0,1)-norm(X'*n,1);
        n = n/norm(n,2);
        n0 = n;
    end
    fprintf('Iter=%d, J=%g\n',k,norm(X'*n0,1))
    disp('Converged.')
    idx0 = find(abs(X'*n) < 1e-3);
    
    if length(idx0) >= 0.3*size(X,2)
        disp('possibly the gobal minimum.')
        B = [B, n];
    else
        disp('possibly a local minimum.')
        B = [];
    end
    X2 = X;
    while (c>=2)
        [U,X2] = projection(X2,n);
        n = DPCP(X2(:,idx0),c-1,e,max_iter);
        if isempty(n)
            disp('terminate.')
            break
        else 
            B = [B, U*n];
        end
        c = c-1;
    end
    
end


