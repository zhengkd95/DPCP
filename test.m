% a simple test
X=[];
B = eye(4);

N=10;
basis = B(:,1:4);
X = [X basis*rand(size(basis,2),N)];

N=100;
basis = B(:,1:3);
X = [X basis*(rand(size(basis,2),N))];

DPCP(X)