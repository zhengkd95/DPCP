% Project the data X onto orthogonal space of vector b
% OUTPUT:
% U: the base vectors of the new subspace
% X2: the dimension reduced data
function [U,X2] = projection(X,b)

d = b'*X;
X2 = X - repmat(d,length(b),1).*repmat(b,1,size(X,2));

[coeff,score,latent,tsquared,explained] = pca(X2','Centered',false); %PCA

idx = find(abs(explained) > 1e-4 );

U = coeff(:,idx);
X2 = score(:,idx)';

end