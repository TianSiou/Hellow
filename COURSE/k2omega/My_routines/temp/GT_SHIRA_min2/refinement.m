function [EV] = refinement(Hessen_H, upper_R, ew, invar_V)

dim_n = length(ew);
EV    = zeros(size(invar_V,1), dim_n);
for ii = 1:dim_n
    mtx       = Hessen_H - ew(ii) * upper_R;
    [U, S, V] = svd(mtx);
    EV(:,ii)  = invar_V * V(:,end);
end