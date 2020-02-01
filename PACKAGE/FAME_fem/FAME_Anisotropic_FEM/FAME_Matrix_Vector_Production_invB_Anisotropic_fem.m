function vec_y = FAME_Matrix_Vector_Production_invB_Anisotropic_fem(vec_x, B)
    if strcmp(B.lssvr, 'mldivide')
            vec_y = B.B\vec_x;
        elseif strcmp(B.lssvr, 'pcg')
%             time_start_pcg = tic;
            vec_y = pcg( B.B, vec_x, 1e-12, 1000);
%             cpu_time_pcg = toc(time_start_pcg);
%             fprintf('Solve linear system Bx=y by pcg using cpu time %3.2f (sec.) with relres = %.2e\n',cpu_time_pcg,relres);
        elseif isfield(B,'P')
            vec_y = B.L\vec_x(B.P);
            vec_y = B.U\vec_y;
            vec_y = vec_y(B.invP);
        else
            vec_y = B.L\vec_x;
            vec_y = B.U\vec_y;
    end
end