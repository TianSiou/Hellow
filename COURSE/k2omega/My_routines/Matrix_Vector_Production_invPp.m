 function vec_y = Matrix_Vector_Production_invPp(vec_x, mode, mtx_palindromic, alpha, fun_Tp, fun_Tps)
    vec_y = fun_Tps(vec_x);
    if strcmp(mode,'transpose') || strcmp(mode,'trans') || strcmp(mode,'transp')
        vec_y = vec_y - (1/alpha)*(mtx_palindromic.Lambda_N2bar*(mtx_palindromic.Lambda_N1bar'*vec_y));
        vec_y = vec_y.*mtx_palindromic.invSigma_pt;
    else
        vec_y = vec_y - (1/alpha)*(mtx_palindromic.Lambda_N1*(mtx_palindromic.Lambda_N2'*vec_y));
        vec_y = vec_y.*mtx_palindromic.invSigma_p;
    end
    vec_y = fun_Tp(vec_y);
end
