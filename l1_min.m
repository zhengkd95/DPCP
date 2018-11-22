% solve the l1 minimum problem:
% n = argmin_{n'*n0=1} ||X'*n||_1
function n = l1_min(X,n0)

settings = sdpsettings('verbose',0);
nk = sdpvar(length(n0),1);
Obj = norm(X'*nk,1);
Constraints = [n0'*nk ==1];

results = optimize(Constraints,Obj,settings);
if results.problem == 0
    n = value(nk);
    n(isnan(n)) = 0;
else
    n = zeros(size(n0));
end

end