function [L_cluster, Perm_vec_all, vec_leng_all, flag] = finding_SiCo_graph_comp(mtx_L,fid)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    s = RandStream('mt19937ar','Seed',0);
    RandStream.setGlobalStream(s);
    
    dim_n          = size(mtx_L,1);
    ew_number      = 3;
    no_restart     = max(100, 2*ew_number);
    stop_tolerance = 1.0e-13 * sqrt(dim_n); %1.0e-12;
%     initial_vec    = rand(dim_n,1);
%     initial_vec    = initial_vec / norm(initial_vec);
    target         = 1.0e-3;
    target_type    = 'CL';
    CorEq          = 'SIRA';
    LSinfo.precond = 'no';
    LSinfo.solver  = 'minres'; %'minres'; %'pcg';
    
    kk       = 1;

    Perm_vec_all = zeros(dim_n,1);
    vec_leng_all = [];
    L_cluster    = mtx_L;
    nn           = dim_n;
    ncluster     = 10;
    flag         = 0;
    flag_s       = 0;
    flag_null    = 1;
    initial_vec  = rand(nn,1);
    initial_vec  = initial_vec / norm(initial_vec);
    
    while ( flag_null == 1 && flag == 0 && flag_s == 0 && kk <= 6 )
    
%         initial_vec    = rand(nn,1);
%         initial_vec    = initial_vec / norm(initial_vec);
    
        [EW, EV, initial_vec, ritz_ew] = NullSp_SEP_Herm_JDSIRA_Driver_v2 (L_cluster, nn, no_restart, ew_number, ...
                stop_tolerance, initial_vec, target, fid, target_type, CorEq, ...
                LSinfo, []);
      
        [~,idx_sort] = sort(abs(EW));
        EW           = EW(idx_sort,1);
        EV           = EV(:,idx_sort);
        
        ncluster     = length(find(abs(EW) < 1.0e-10));
        
        
        if ( ncluster == 1 )
            err = sign(max(EV(:,1))) * EV(:,1) - ones(nn,1) / sqrt(nn); 
            if (max(abs(err)) > stop_tolerance*1.0e3)
                flag_null   = 1; 
            else
                flag_null   = 0;
            end
        elseif ( ncluster > 1 ) 
            flag_null   = 1;
        else
            flag_null   = 0;
        end
  
        %[flag_null ncluster min(err) max(err)]
                
        if ( flag_null == 1 )
            
            %fprintf('ncluster = %3.0f, ncluster_estimate + ncluster_1 = %3.0f, ncv = %3.0f \n', ncluster, ncluster_estimate + ncluster_1, ncv);
        
            tic;
            [Perm_vec, L_cluster_tmp, vec_leng, flag] = one_null_cluster(L_cluster, EV(:,1)); 
            
%             [Perm_vec, L_cluster_tmp, vec_leng, flag] = kmeans_cluster(L_cluster, real(EV(:,1:ncluster)), ncluster);
            fprintf(fid,'CPU time for clustering is %11.4e and number of clusters is %7.0f. \n\n', toc, length(vec_leng));
            
            if ( ncluster == 1 || ritz_ew(end,1) > 1.0e-3 )
                initial_vec = rand(vec_leng(1),1);
                initial_vec = initial_vec / norm(initial_vec);
            else
                initial_vec      = initial_vec(Perm_vec,:);
                initial_vec      = initial_vec(1:vec_leng(1),:);
                [~, initial_vec] = Gram_Schmidt(size(initial_vec,2), initial_vec);
            end
            
            stop_tolerance = 1.0e-10;
%             target         = -1.0e-3;
%             LSinfo.solver  = 'pcg';
        
            if ( kk == 1 )
                Perm_vec_all = Perm_vec;
            else
                tmp                  = Perm_vec_all(1:nn,1);
                Perm_vec_all(1:nn,1) = tmp(Perm_vec);
            end 
            
            if ( flag == 0 )
                L_s        = L_cluster_tmp(vec_leng(1,1)+1:end,vec_leng(1,1)+1:end);
                ns         = size(L_s,1);
                if ( ns < 500 ) %5000 )
                    tic;
%                     [EVs, EWs] = eig(full(L_s));
%                     EWs        = diag(EWs);
                    EWs        = eig(full(L_s));
                    ncluster_s = length(find(abs(EWs) < 1.0e-10));
%                     [Perm_vec_s, L_cluster_s, vec_leng_s, flag_s] = kmeans_cluster(L_s, real(EVs(:,1:ncluster_s)), ncluster_s);
                    fprintf(fid,'%2.0f-th, number of zero eigenvalues of the 2nd cluster matrix is %3.0f, CPU time = %11.4e \n', kk, ncluster_s, toc);
                    fprintf('%2.0f-th, number of zero eigenvalues of the 2nd cluster matrix is %3.0f \n', kk, ncluster_s);
                    
                    vec_leng_s = vec_leng(2:end,1);
                    
                    if ( ncluster_s == 0 )
                        flag_s = 1;
                    else
                        flag_s = 0;
                    end
                    
%                     tmp                                = Perm_vec_all(vec_leng(1,1)+1:nn,1); 
%                     Perm_vec_all(vec_leng(1,1)+1:nn,1) = tmp(Perm_vec_s,1);
                
                else
                    
                    fprintf(fid,'Matrix dimension of 2nd cluster matrix is %7.0f \n\n',ns);
%                     [L_cluster_s, Perm_vec_s, vec_leng_s, flag_s] = finding_SiCo_graph_eigs(L_s,fid);
                    idx_end = vec_leng(1);
                    zero_ew = 0;
                    tic;
                    for ii = 2:length(vec_leng)
                        idx_begin  = idx_end + 1;
                        idx_end    = idx_end + vec_leng(ii);
                        EWs        = eig(full(L_cluster_tmp(idx_begin:idx_end,idx_begin:idx_end))); 
                        ncluster_s = length(find(abs(EWs) < 1.0e-10));
                        
                        if ( ncluster_s == 0 )
                            flag_s = 1;
                            fprintf(fid,'Error in checking number of zero eigenvalues of the 2nd cluster matrix \n');
                            fprintf('Error in checking number of zero eigenvalues of the 2nd cluster matrix \n');
                            save testing_data L_cluster_tmp mtx_L Perm_vec_all vec_leng Perm_vec
                        else
                            zero_ew = zero_ew + ncluster_s;
                            flag_s  = 0;
                        end
                    end
                    
                    vec_leng_s = vec_leng(2:end,1);
                     
                    fprintf(fid,'%2.0f-th, number of zero eigenvalues of the 2nd cluster matrix is %3.0f with CPU time %11.4e \n', kk, zero_ew, toc);
                    fprintf('%2.0f-th, number of zero eigenvalues of the 2nd cluster matrix is %3.0f \n', kk, zero_ew);
                    
                end
        
%                 tmp                                = Perm_vec_all(vec_leng(1,1)+1:nn,1); 
%                 Perm_vec_all(vec_leng(1,1)+1:nn,1) = tmp(Perm_vec_s,1);
                
                if ( flag_s == 0 )
                    nn           = vec_leng(1,1);
                    L_cluster    = L_cluster_tmp(1:nn,1:nn);
            
                    vec_leng_all = [ vec_leng_s; vec_leng_all];
                end
            
            end
            kk        = kk + 1; 
    
        end
    
    end
    
    if ( ncluster == 1 && kk > 1 ) % success to cluster
        vec_leng_all = [ vec_leng(1,1); vec_leng_all];
        flag         = 0;
    elseif ( ncluster == 0 && norm(L_cluster * ones(size(L_cluster,1),1)) < 1.0e-12 ) % success to cluster
        vec_leng_all = [ vec_leng(1,1); vec_leng_all];
        flag         = 0;
    else % fail to cluster
        flag = 1;
    end

    if ( ncluster == 1 && kk == 1 && flag_null == 0 )
        L_cluster    = mtx_L;
        flag         = 0;
        vec_leng_all = dim_n;
        
        fprintf(fid, 'number of zero eigenvalues of the Laplacian matrix is 1 \n');
        fprintf('number of zero eigenvalues of the Laplacian matrix is 1 \n'); 
    
    else
        L_cluster = mtx_L(Perm_vec_all, Perm_vec_all);
        
        fprintf(fid, 'number of zero eigenvalues of the Laplacian matrix is %3.0f \n',length(vec_leng_all));
        fprintf('number of zero eigenvalues of the Laplacian matrix is %3.0f \n',length(vec_leng_all));
        fprintf(fid, 'dimesion of 1st cluster is %8.0f \n', vec_leng_all(1,1));
        fprintf('dimesion of 1st cluster is %8.0f \n', vec_leng_all(1,1));
    end

end

