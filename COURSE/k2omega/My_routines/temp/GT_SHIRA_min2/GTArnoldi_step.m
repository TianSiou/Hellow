function [Y, Z, Htilde, Rtilde] = GTArnoldi_step(n, A, Q, solve_LS_PQEP, old_Y, old_Z, H, R, shift, loopstart, loopend)
    Y                   = zeros(2*n, loopend+1);
    Z                   = zeros(2*n, loopend);
    Y(:, 1:loopstart)   = old_Y;
    Z(:, 1:loopstart-1) = old_Z;
    rhat                = zeros(loopend,1);
    
    Z_Htran                   = zeros(loopend, 2*n);
    Y_Htran                   = zeros(loopend+1, 2*n);
    Y_Htran(1:loopstart, :)   = old_Y';
    Z_Htran(1:loopstart-1, :) = old_Z';
    Y_conj                    = zeros(2*n, loopend+1);
    Z_conj                    = zeros(2*n, loopend);
    Y_conj(:, 1:loopstart)    = conj(old_Y);
    Z_conj(:, 1:loopstart-1)  = conj(old_Z);
    
    Rtilde                               = zeros(loopend);
    Rtilde(1:loopstart-1, 1:loopstart-1) = R;
    Htilde                               = zeros(loopend+1,loopend);
    Htilde(1:loopstart,   1:loopstart-1) = H;
    
%
    % Steps of Arnoldi process from loopstart to loopend.
    %
    for j = loopstart:loopend
    %
    %  Solve linear system
    %       Ntilde*Z(:,j) = Y(:,j),
    %  where 
    %       Ntilde = ( M - shift L ) J ( M^T - shift L^T ) J^T.
    % 
       [Z(:,j)]     = GT_linear_system(n, A, Q, solve_LS_PQEP, Y(:,j), shift);
       %Z_Htran(j,:) = Z(:,j)';
       
%        M = [A, zeros(n, n); Q, -eye(n)];
%        L = [zeros(n, n), eye(n); A.', zeros(n, n)];
%        J = [zeros(n, n), eye(n); -eye(n), zeros(n, n)];
%        N = (M - shift*L)*J*((M.')-shift*(L.'))*(J.');
%        
%        Z(:, j) = N\Y(:, j);
%        y = mtx_N_prod_vec( A, Q, Z(:,j), shift );
%        norm(Y(:, j)-y)
       
    %
    % full reorthogonal Z(:,j) to previous Z(:,1:j-1)
    %   
    
        for ii = 1:j-1
            rhat(ii,1) = Z_Htran(ii,:) * Z(:,j);
            Z(:,j)     = Z(:,j) - rhat(ii,1) * Z(:,ii);
        end
        
%     %
%     % reorthogonalize against Z
%     %
%     
%         for ii = 1:j-1
%            cor        = Z(:,ii)' * Z(:,j);
%            Z(:,j)     = Z(:,j) - cor * Z(:,ii);
%            rhat(ii,1) = rhat(ii,1) + cor;
%         end
%         
    %
    % full resymplectic Z(:,j) to Y(:.1:j)
    % 
    
        for ii = 1:j
            tmp_s  = [ Y(n+1:n*2,ii).'  -Y(1:n,ii).' ] * Z(:,j);
            Z(:,j) = Z(:,j) - tmp_s * [ Y_conj(n+1:n*2,ii);  -Y_conj(1:n,ii) ];
        end
        
        Rtilde(j,j)     = 1 / norm(Z(:,j));
        Z(:,j)          = Rtilde(j,j) * Z(:,j);
        Z_Htran(j,:)    = Z(:,j)';
        Z_conj(:,j)     = conj(Z(:,j));
        rhat(1:j-1,1)   = Rtilde(j,j) * rhat(1:j-1,1);
        Rtilde(1:j-1,j) = - Rtilde(1:j-1,1:j-1) * rhat(1:j-1,1);
        
    %
    %   Compute K Z(:,j)
    %
    
%         zj = Z(:, j);
%         K = [A, zeros(n, n); zeros(n, n), A.'];
%         y = K*Z(:, j);
%         y = y*shift;
        Y(:, j+1)  =  mtx_K_prod_vec( A, Z(:, j), shift );
%         tmp            = A*Z(1:n, j);
%         Y(1 :n, j+1)   = shift * tmp;
%         tmp            = A.'*Z((n+1):(2*n), j);
%         Y((n+1):(2*n),j+1) = shift * tmp;
%         norm(y - Y(:, j+1))
        
    %
    % full reorthogonal Y(:,j+1) to previous Y(:,1:j)
    %
    
        for ii = 1:j
            Htilde(ii,j) = Y_Htran(ii,:) * Y(:,j+1);
            Y(:,j+1)     = Y(:,j+1) - Htilde(ii,j) * Y(:,ii);
        end
%         norm(Y(:, j+1).'*Y(:, 1:j))
        
%     %
%     % reorthogonalize against Y
%     %
%         for ii = 1:j
%            cor          = Y(:,ii)' * Y(:,j+1);
%            Y(:,j+1)     = Y(:,j+1) - cor * Y(:,ii);
%            Htilde(ii,j) = Htilde(ii,j) + cor;
%         end
% %         norm(Y(:, 1:j).'*Y(:, j+1))
        
    %
    % resymplectic Y(:,j+1) to Z(:.1:j)
    %
    
        for ii = 1:j
            tmp_t    = [ Z(n+1:2*n,ii).'  -Z(1:n,ii).' ] * Y(:,j+1);
            Y(:,j+1) = Y(:,j+1) - tmp_t * [ Z_conj(n+1:2*n,ii);  -Z_conj(1:n,ii) ];
        end
%         norm(Y(:, j+1).'*J*conj(Z(:, 1:j)))
        
        Htilde(j+1,j)  = norm(Y(:,j+1));
        Y(:,j+1)       = Y(:,j+1) / Htilde(j+1,j);
        Y_Htran(j+1,:) = Y(:,j+1)';
        Y_conj(:,j+1)  = conj(Y(:,j+1));
%         norm(Y(:, 1:j+1).'*Y(:, 1:j+1) - eye(j+1))
    end
end