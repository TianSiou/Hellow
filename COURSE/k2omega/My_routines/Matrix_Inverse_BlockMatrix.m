function invD_fun = Matrix_Inverse_BlockMatrix(D)
% Construct inverse of a block matrix
%           [ diag(D.d11), diag(D.d12), diag(D.d13) ]
%       D = [ diag(D.d21), diag(D.d22), diag(D.d33) ]
%           [ diag(D.d31), diag(D.d32), diag(D.d33) ]
%     N = size(D,1)/3;
%     D11 = diag(D(1:N,    1:  N));
%     D12 = diag(D(1:N,  N+1:2*N));
%     D13 = diag(D(1:N,2*N+1:3*N));
%     
%     D21 = diag(D(N+1:2*N,    1:  N));
%     D22 = diag(D(N+1:2*N,  N+1:2*N));
%     D23 = diag(D(N+1:2*N,2*N+1:3*N));
%     
%     D31 = diag(D(2*N+1:3*N,    1:  N));
%     D32 = diag(D(2*N+1:3*N,  N+1:2*N));
%     D33 = diag(D(2*N+1:3*N,2*N+1:3*N));
    
    
%     N = length(D11);
    
%     det    = -D11.*D22.*D33 + D11.*D23.*D32 + D12.*D21.*D33 - D12.*D31.*D23 - D21.*D13.*D32 + D13.*D22.*D31;
    det    = -D.d11.*D.d22.*D.d33 + D.d11.*D.d23.*D.d32 + D.d12.*D.d21.*D.d33 - D.d12.*D.d31.*D.d23 - D.d21.*D.d13.*D.d32 + D.d13.*D.d22.*D.d31;
    invdet = @(vec_x) [ vec_x(1:end/3)./det; vec_x(end/3+1:2*end/3)./det; vec_x(2*end/3+1:end)./det ];

%     invD_11 = (-D22.*D33 + D23.*D32)./det;
%     invD_22 = (-D11.*D33 + D13.*D31)./det;
%     invD_33 = (-D11.*D22 + D12.*D21)./det;
%     invD_12 = ( D12.*D33 - D13.*D32)./det;
%     invD_21 = ( D21.*D33 - D31.*D23)./det;
%     invD_13 = (-D12.*D23 + D13.*D22)./det;
%     invD_31 = (-D21.*D32 + D22.*D31)./det;
%     invD_23 = ( D11.*D23 - D21.*D13)./det;
%     invD_32 = ( D11.*D32 - D12.*D31)./det;
    temp.temp_invD_11 = (-D.d22.*D.d33 + D.d23.*D.d32);
    temp.temp_invD_22 = (-D.d11.*D.d33 + D.d13.*D.d31);
    temp.temp_invD_33 = (-D.d11.*D.d22 + D.d12.*D.d21);
    temp.temp_invD_12 = ( D.d12.*D.d33 - D.d13.*D.d32);
    temp.temp_invD_21 = ( D.d21.*D.d33 - D.d31.*D.d23);
    temp.temp_invD_13 = (-D.d12.*D.d23 + D.d13.*D.d22);
    temp.temp_invD_31 = (-D.d21.*D.d32 + D.d22.*D.d31);
    temp.temp_invD_23 = ( D.d11.*D.d23 - D.d21.*D.d13);
    temp.temp_invD_32 = ( D.d11.*D.d32 - D.d12.*D.d31);
    
    temp_invD = @(vec_x) [ temp.temp_invD_11.*vec_x(1:end/3) + temp.temp_invD_12.*vec_x(end/3+1:2*end/3) + temp.temp_invD_13.*vec_x(2*end/3+1:end);
                           temp.temp_invD_21.*vec_x(1:end/3) + temp.temp_invD_22.*vec_x(end/3+1:2*end/3) + temp.temp_invD_23.*vec_x(2*end/3+1:end);
                           temp.temp_invD_31.*vec_x(1:end/3) + temp.temp_invD_32.*vec_x(end/3+1:2*end/3) + temp.temp_invD_33.*vec_x(2*end/3+1:end) ];
    invD_fun = @(vec_x) temp_invD( invdet(vec_x) );

%     D = sparse(    1:3*N,     1:3*N, [D11;D22;D33], 3*N,3*N) + ...
%         sparse(  N+1:3*N,     1:2*N,     [D21;D32], 3*N,3*N) + ...
%         sparse(2*N+1:3*N,     1:  N,           D31, 3*N,3*N) + ...
%         sparse(    1:2*N,   N+1:3*N,     [D12;D23], 3*N,3*N) + ...
%         sparse(    1:  N, 2*N+1:3*N,           D13, 3*N,3*N);
%     invD.mtx = sparse(    1:3*N,     1:3*N, [invD_11;invD_22;invD_33], 3*N,3*N) + ...
%                sparse(  N+1:3*N,     1:2*N,         [invD_21;invD_32], 3*N,3*N) + ...
%                sparse(2*N+1:3*N,     1:  N,                   invD_31, 3*N,3*N) + ...
%                sparse(    1:2*N,   N+1:3*N,         [invD_12;invD_23], 3*N,3*N) + ...
%                sparse(    1:  N, 2*N+1:3*N,                   invD_13, 3*N,3*N); 
end

function fun = fun_D(vec_x, temp)
    fun = @(vec_x) [ temp.temp_invD_11.*vec_x(1:end/3) + temp.temp_invD_12.*vec_x(end/3+1:2*end/3) + temp.temp_invD_13.*vec_x(2*end/3+1:end);
                     temp.temp_invD_21.*vec_x(1:end/3) + temp.temp_invD_22.*vec_x(end/3+1:2*end/3) + temp.temp_invD_23.*vec_x(2*end/3+1:end);
                     temp.temp_invD_31.*vec_x(1:end/3) + temp.temp_invD_32.*vec_x(end/3+1:2*end/3) + temp.temp_invD_33.*vec_x(2*end/3+1:end) ];
end

% function invD = Matrix_Inverse_BlockMatrix(D)
% % Construct inverse of a block matrix
% %           [ diag(D11), diag(D12), diag(D13) ]
% %       D = [ diag(D21), diag(D22), diag(D23) ]
% %           [ diag(D31), diag(D32), diag(D33) ]
%     N = size(D,1)/3;
%     D11 = diag(D(1:N,    1:  N));
%     D12 = diag(D(1:N,  N+1:2*N));
%     D13 = diag(D(1:N,2*N+1:3*N));
%     
%     D21 = diag(D(N+1:2*N,    1:  N));
%     D22 = diag(D(N+1:2*N,  N+1:2*N));
%     D23 = diag(D(N+1:2*N,2*N+1:3*N));
%     
%     D31 = diag(D(2*N+1:3*N,    1:  N));
%     D32 = diag(D(2*N+1:3*N,  N+1:2*N));
%     D33 = diag(D(2*N+1:3*N,2*N+1:3*N));
%     
%     
%     N = length(D11);
%     
%     det = -D11.*D22.*D33 + D11.*D23.*D32 + D12.*D21.*D33 - D12.*D31.*D23 - D21.*D13.*D32 + D13.*D22.*D31;
% 
%     invD_11 = (-D22.*D33 + D23.*D32)./det;
%     invD_22 = (-D11.*D33 + D13.*D31)./det;
%     invD_33 = (-D11.*D22 + D12.*D21)./det;
%     invD_12 = ( D12.*D33 - D13.*D32)./det;
%     invD_21 = ( D21.*D33 - D31.*D23)./det;
%     invD_13 = (-D12.*D23 + D13.*D22)./det;
%     invD_31 = (-D21.*D32 + D22.*D31)./det;
%     invD_23 = ( D11.*D23 - D21.*D13)./det;
%     invD_32 = ( D11.*D32 - D12.*D31)./det;
% 
% %     D = sparse(    1:3*N,     1:3*N, [D11;D22;D33], 3*N,3*N) + ...
% %         sparse(  N+1:3*N,     1:2*N,     [D21;D32], 3*N,3*N) + ...
% %         sparse(2*N+1:3*N,     1:  N,           D31, 3*N,3*N) + ...
% %         sparse(    1:2*N,   N+1:3*N,     [D12;D23], 3*N,3*N) + ...
% %         sparse(    1:  N, 2*N+1:3*N,           D13, 3*N,3*N);
%     invD.mtx = sparse(    1:3*N,     1:3*N, [invD_11;invD_22;invD_33], 3*N,3*N) + ...
%                sparse(  N+1:3*N,     1:2*N,         [invD_21;invD_32], 3*N,3*N) + ...
%                sparse(2*N+1:3*N,     1:  N,                   invD_31, 3*N,3*N) + ...
%                sparse(    1:2*N,   N+1:3*N,         [invD_12;invD_23], 3*N,3*N) + ...
%                sparse(    1:  N, 2*N+1:3*N,                   invD_13, 3*N,3*N); 
% end