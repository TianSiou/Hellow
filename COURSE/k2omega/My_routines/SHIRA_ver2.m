function varargout = SHIRA_ver2( varargin )
%SHIRA  Find desired eigenpairs of skew-Hamiltonian matrix A.
%   D = SHIRA(A) returns a vector of 6 largest magnitude eigenvalues.
%   A must be skew-Hamiltonian.
%
%   [V,D] = SHIRA(A) returns a diagonal matrix D of 6 largest magnitude
%   eigenvalues and a matrix V whose columns are the corresponding eigenvectors.
%
%   [V,D,R] = SHIRA(A) also returns a column vector of corresponding residuals. 
%
%   SHIRA(A,k) return the k largest magnitude eigenvalues.
%
%   SHIRA(A,k,tau) return k eigenvalues. If tau is:
%      'lm' or 'sm' - Largest or Smallest Magnitude
%   If tau is a real or complex scalar including 0, SHIRA finds these
%   eigenvalues which closest to tau.
    global check_add_eigen_wanted
    check_add_eigen_wanted = 0;
    [A , eigen_wanted, sigma, option ] = input_check(nargin,nargout,varargin);
   
       %% Set parameters of Arnoldi step
    global n

    start_vec    = option.v0;
    remax        = 20000;                 % max restart number
%     eigen_wanted = ceil(eigen_wanted/2);  % sufficient to determine half
    
    
       %% SHIRA process
    iter   = 0;   % total iterations
    nre    = 0;   % restart number
    m      = max([3*eigen_wanted , 30 ]); % size of subspace generated by arnoldi step 
    k      = 0;  % size of subspace after restart
    U      = zeros(n,m);
    H      = zeros(m,m);
    U(:,1) = start_vec; 
    while nre < remax
        nre = nre + 1;
        
        [ U, H ] = isotropic_Arnoldi_step( A , sigma, k , m , U , H, option );
        iter = m - k + iter;
        
        [ U, H, hessenberg_residual, k, eigen_wanted ] = implicitly_restart_step(U, H, eigen_wanted); % In this version, k = eigen_wanted
hessenberg_residual
        if  hessenberg_residual < (1e-14)*( abs(H(k+1,k+1)) + abs(H(k,k)) ) 
            U = U(:,1:k);
            H = H(1:k,1:k);
            
            [eig_vec, eig_val] = eig(H);
            % recovery eigenpairs
            if strcmp(sigma,'sm')
                D = diag(1./diag(eig_val));
            elseif isnumeric(sigma)==1
                D = diag(1./diag(eig_val) + sigma);
            else
                D       = eig_val;
            end
            eig_vec = U*eig_vec;
            V       = eig_vec;
   
            varargout = output_check( V, D, U, H, nargout );
            return;
        end
        
%         if nre > 1 && hessenberg_residual > hessenberg_residual_old  && hessenberg_residual_old < 1e-8
%             U(:,2:m) = 0;
%             H        = zeros(m,m);
%             U(:,1)   = U(:,1) + 1e-2;
%             U(:,1)   = U(:,1)/norm(U(:,1));
%             k = 0;
%             check_add_eigen_wanted = 0;
%         end
%         hessenberg_residual_old = hessenberg_residual;
        
    end
end

%% Input check
function [ A, eigen_wanted, sigma, option ] = input_check(Nargin,Nargout,Varargin)
    global handle_mode n
    
    eigen_wanted  = 6;
    sigma         = 'lm';
    option = struct(            'v0', [],...
                     'linear_solver', 'mldivide',...
                            'isreal', 0);
    switch Nargin
        case 0 
            error ('Not enough input!');
        case 1
            A             = Varargin{1};
            if isvector(A) 
                n = Varargin{2};
            else
                n = length(A);
            end
        case 2
            A             = Varargin{1};
            if isvector(A) 
                n = Varargin{2};
            else
                n = length(A);
                eigen_wanted  = Varargin{2};
            end
        case 3
            A             = Varargin{1};
            if isvector(A) 
                n             = Varargin{2};
                eigen_wanted  = Varargin{3};
            else
                n             = length(A);
                eigen_wanted  = Varargin{2};
                sigma         = Varargin{3};
            end
        case 4
            A             = Varargin{1};
            if isvector(A) 
                n = Varargin{2};
                eigen_wanted  = Varargin{3};
                sigma         = Varargin{4};
            else
                n = length(A);
                eigen_wanted  = Varargin{2};
                sigma         = Varargin{3};
                option        = Varargin{4};
            end
        case 5
            A             = Varargin{1};
            if isvector(A) 
                n = Varargin{2};
                eigen_wanted  = Varargin{3};
                sigma         = Varargin{4};
                option        = Varargin{5};
            else
                error('Too many input argument')
            end
        otherwise
            error('Too many input argument')
    end
    handle_mode = isvector(A); % handle string:1, matrix data:0
    if isempty(option.v0)
        option.v0 = rand(n,1);
        option.v0 = option.v0/norm(option.v0);
    end
       
    eigen_wanted = floor(eigen_wanted);
    if isnumeric(sigma)==1 && length(sigma)~=1
        error('The nearby target tau must be a real or complex number');
    elseif isnumeric(sigma)==0
        switch sigma
            case 'lm'
                sigma = 'lm';
            case 'sm'
                sigma = 'sm';
            otherwise
                error('Invalid input of tau');
        end
    end
    if Nargout > 5
        error('Too many output argument.');
    end
end
%% Input check
function Varargout = output_check( EV, EW, U, H, Nargout )
    switch Nargout
        case 0 || 1
            Varargout{1} = diag(EW);
        case 2
            Varargout{1} = EV;
            Varargout{2} = EW;
        case 3
            Varargout{1} = EV;
            Varargout{2} = EW;
            Varargout{3} = U;
        case 4
            Varargout{1} = EV;
            Varargout{2} = EW;
            Varargout{3} = U;
            Varargout{4} = H;
        case 5
            Varargout{1} = EV;
            Varargout{2} = EW;
            Varargout{3} = U;
            Varargout{4} = H;
            Varargout{5} = iter;
    end
end
%% isotropic Arnoldi step
function [ U, H ] = isotropic_Arnoldi_step( A, sigma, k, m, U, H, option )
    % span the isotropic Arnoldi decomposition from k dimension to m dimension
    global handle_mode n

    for j = k+1:m
        if handle_mode == 0
            switch sigma
                case 'lm'
                    temp = A*U(:,j);
                case 'sm'
                    temp = Mtx_vec_product_A_shift_I( A, 0, U(:,j), option );
                otherwise
                    temp = Mtx_vec_product_A_shift_I( A, sigma, U(:,j), option );
            end
        else
            temp = A(U(:,j));
        end
        
        UU = [U(1:n/2,1:j)  -U(n/2+1:n,1:j) ; U(n/2+1:n,1:j)  U(1:n/2,1:j)];  % UU = [U,-JU];
        [u, h, beta] = gsreorthog(UU , temp);  % orthogonal temp against [U,-JU]
%         [u, h, beta] = gsreorthog(full([U(1:n/2,1:j)  -U(n/2+1:n,1:j) ; U(n/2+1:n,1:j)  U(1:n/2,1:j)]) , temp);  % orthogonal temp against [U,-JU]
        
        H(1:j,j)     = h(1:j);
        H(j+1,j)     = beta;
        U(:,j+1)     = u;
    end
end
function y = Mtx_vec_product_A_shift_I( A, sigma, x, option )
global n
    switch option.linear_solver
        case 'mldivide'
            y = (A-sigma*eye(n))\x;    
    end

end
%% reorthogonalize process
function [p, r, tho] = gsreorthog(Q, x)
% [p, r, tho] = gsreorthog(Q, x)
%      x = Q*r + tho*p 
% Orthogonalize x against Q with reorthogonalization
    [m ,n] = size(Q);
%     r = zeros(n,1);
    r = zeros(n,1);
    tho = norm(x);
    if tho == 0
%         p = zeros(n,1);
        p = zeros(n,1);
        return
    end
    mu = tho;
    x_bar = x;
    
    % find the index of the column with minimal 1-norm
    n1 = norm(Q(:,1));
    k = 1;
    for i=2:n
        n2 = norm(Q(:,i));
        if(n1>n2)
            n1 = n2; 
            k = i;
        end
    end

%     ei = zeros(m,1);
    ei = zeros(m,1);
    ei(k) = 1;
    temp = Q';
    while(1) 
        s = temp*x_bar;
%         s = Q'*x_bar;
        r = r + s;
        
        % x_bar is the residual
        x_bar = x_bar - Q*s;
        tau = norm(x_bar);
        if(tau/tho >=1/2)
            break;
        end
       
        if(tau > 0.1*mu*eps)
            tho = tau;
        else
            tho = 0.1*tho*eps;
            mu = tho;
            x_bar = ei;
            display('reorthogonal')
        end
    end
    tho = norm(x_bar);
    p = x_bar/tho;
end
%% Implicitly restarted process
function [ U , H , hessenberg_residual , k, eigen_wanted] = implicitly_restart_step(U , H , eigen_wanted )
    % reduce the m dimensional isotropic Arnoldi decomposition to the k dimensional.
    m  =  size(U,2) - 1;
        
    beta =  H(m+1,m);
    U_m  =  U(:,1:m);       
    H_m  =  H(1:m,1:m);
    
    [shift_eig, eigen_wanted] = choose_shift( H_m, beta, eigen_wanted);
    
    [U_m , H_m , b] = QR_step(U_m, H_m, shift_eig);
    
    k      = m - length(shift_eig);
    w      = H_m(k+1,k)*U_m(:,k+1) + beta*b(k)*U(:,m+1);
    beta   = norm(w);
    u      = w/beta;
    
    U(:,1:k+1) = [U_m(:,1:k) u];
    H(1:k,1:k) = H_m(1:k,1:k);
    H(k+1,k)   = beta;
    
    hessenberg_residual = abs(beta);
end
%% choose unwanted ritz values
function [ shift_eig, eigen_wanted ] = choose_shift( H_m, beta, eigen_wanted)  
    % choose unwanted ritz values from H_m
    global check_add_eigen_wanted
    
    m = length(H_m);
    
    [Y , ritz_val] = eig(H_m);
% Res_ritz = abs(beta*Y(end,:))
    ritz_val = diag(ritz_val);
    
    % Choose strategy 1 (in terms of residual)
%     residual = beta*Y(m,:);  % compute residual of approximate ritz pairs  
%     [~ , S]  = sort(residual,'ascend'); 
    % Choose strategy 2 (in terms of magnitude)
    magnitude = abs(ritz_val);
    [~ , S]  = sort(magnitude,'descend'); 
    
    ritz_val = ritz_val(S);
    
    if abs( ritz_val(eigen_wanted) - conj(ritz_val(eigen_wanted+1)) ) < 1e-8  &&  check_add_eigen_wanted == 0
        eigen_wanted = eigen_wanted + 1;
        check_add_eigen_wanted = 1;
    elseif abs( ritz_val(eigen_wanted) - conj(ritz_val(eigen_wanted+1)) ) < 1e-8  &&  check_add_eigen_wanted == 1
        eigen_wanted = eigen_wanted - 1;
        check_add_eigen_wanted = 0;
    end

    shift_eig = ritz_val(eigen_wanted+1 : m);
end
%% Do the QR step to shift unwanted ritz value
function [U_m , H_m , b] = QR_step(U_m, H_m, shift_eig)
    % 
    m = length(H_m);
    e_m  = [zeros(m-1,1) ; 1 ];
    b    = e_m;
    if length(shift_eig) == 1
        [U_m , H_m , b] = single_shift( U_m , H_m , b , shift_eig(1));
    else
        for i = 1:length(shift_eig)
            if i == 1
                if  abs( real(shift_eig(i)-shift_eig(i+1)) ) < 1e-8 && abs( imag(shift_eig(i)+shift_eig(i+1)) ) < 1e-8
                    [U_m , H_m , b] = double_shift( U_m , H_m , b , shift_eig(i));
                else
                    [U_m , H_m , b] = single_shift( U_m , H_m , b , shift_eig(i));
                end
            elseif i ~= length(shift_eig)
                if abs( real(shift_eig(i)-shift_eig(i+1)) ) < 1e-8 && abs( imag(shift_eig(i)+shift_eig(i+1)) ) < 1e-8
                    [U_m , H_m , b] = double_shift( U_m , H_m , b , shift_eig(i));
                elseif abs( real(shift_eig(i)-shift_eig(i-1)) ) < 1e-8 && abs( imag(shift_eig(i)+shift_eig(i-1)) ) < 1e-8
                    continue
                else
                    [U_m , H_m , b] = single_shift( U_m , H_m , b , shift_eig(i));
                end
            else
                if abs( real(shift_eig(i)-shift_eig(i-1)) ) < 1e-8 && abs( imag(shift_eig(i)+shift_eig(i-1)) ) < 1e-8
                    continue
                else
                    [U_m , H_m , b] = single_shift( U_m , H_m , b , shift_eig(i));
                end
            end
        end
    end
end

%% Double shift
function [U , H , b] = double_shift( U , H , b , shift_eig)
    s = 2*real(shift_eig);
    t = real(shift_eig)^2 + imag(shift_eig)^2;
    M = H^2 - s*H + t*eye(size(H));
    
    [Z,~] = qr(M);
    H     = Z'*H*Z;
    H     = triu(H) + diag(diag(H,-1),-1);
    U     = U*Z;
    b     = Z'*b;
end

%% Single shift
function [U , H , b] = single_shift( U , H , b , shift_eig)
    [Q,~] = qr(H - shift_eig*eye(size(H)));
    H     = Q'*H*Q;
    H     = triu(H) + diag(diag(H,-1),-1);
    U     = U*Q;
    b     = Q'*b;
end