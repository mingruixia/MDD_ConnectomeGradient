function [emb,res] = x_compute_diffusion_map(L,alpha,n_components)

if nargin < 2
    alpha = 0.5;
end

ndim = length(L);
L_alpha = L;

d = sum(L_alpha);
d_alpha = d .^ -alpha;

L_alpha = L_alpha.*(d_alpha'*d_alpha);

d_alpha = sum(L_alpha) .^ -1;
L_alpha = repmat(d_alpha',1,ndim) .* L_alpha;
M = L_alpha;

[vectors,lambdas] = eig(M);
lambdas = lambdas(lambdas~=0)';
clear M

[lambdas,lambda_idx] = sort(lambdas, 'descend');
vectors = vectors(:,lambda_idx);

n = max(2, floor(sqrt(ndim))); 
lambdas = lambdas(1:n);
vectors = vectors(:,1:n);


psi = vectors ./ repmat(vectors(:,1),1,size(vectors,2));

diffusion_times = exp(1 - log(1 - lambdas(2:end))./log(lambdas(2:end)));
lambdas = lambdas(2:end) ./ (1 - lambdas(2:end));
lambda_ratio = lambdas ./ lambdas(1);
threshold = max(0.05, lambda_ratio(end));
n_components_auto = min([sum(lambda_ratio>threshold)-1, ndim]);

if nargin < 3
    n_components = n_components_auto;    
end

emb = psi(:,2:n_components + 1) .* repmat(lambdas(1:n_components),ndim,1);
res.lambdas = lambdas;
res.vectors = vectors;
res.n_components = n_components;
res.diffusion_time = diffusion_times;
