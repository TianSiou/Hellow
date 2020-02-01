function [ X, Y ] = J_orthogonal( X, Y, n )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

X(:,1) = X(:,1) / norm(X(:,1));
Y(:,1) = Y(:,1) / (X(:,1)' * [ Y(n+1:2*n,1); -Y(1:n,1) ]);

alpha  = Y(:,1)'*[ X(n+1:2*n,2); -X(1:n,2)];
beta   = X(:,1)'*[ -X(n+1:2*n,2); X(1:n,2)];
X(:,2) = X(:,2) + alpha * X(:,1) + beta * Y(:,1);

alpha  = Y(:,1)'*[ Y(n+1:2*n,2); -Y(1:n,2)];
beta   = X(:,1)'*[ -Y(n+1:2*n,2); Y(1:n,2)];
Y(:,2) = Y(:,2) + alpha * X(:,1) + beta * Y(:,1);

% tt    = X(:,1)'*Y(:,1);
% mtx_A = [ X(:,1)'*X(:,1), tt; tt, Y(:,1)'*Y(:,1)];
% mtx_B = [ X(:,1)'*[ -X(n+1:2*n,2); X(1:n,2)], X(:,1)'*[ -Y(n+1:2*n,2); Y(1:n,2)]; Y(:,1)'*[ -X(n+1:2*n,2); X(1:n,2)], Y(:,1)'*[ -Y(n+1:2*n,2); Y(1:n,2)] ];
% mtx_C = mtx_A \ mtx_B;
% 
% X(:,2) = X(:,2) + mtx_C(1,1) * [ -X(n+1:2*n,1); X(1:n,1) ] + mtx_C(2,1) * [ -Y(n+1:2*n,1); Y(1:n,1) ];
% Y(:,2) = Y(:,2) + mtx_C(1,2) * [ -X(n+1:2*n,1); X(1:n,1) ] + mtx_C(2,2) * [ -Y(n+1:2*n,1); Y(1:n,1) ];

X(:,2) = X(:,2) / norm(X(:,2));
Y(:,2) = Y(:,2) / (X(:,2)'*[ Y(n+1:2*n,2); -Y(1:n,2)]);

end

