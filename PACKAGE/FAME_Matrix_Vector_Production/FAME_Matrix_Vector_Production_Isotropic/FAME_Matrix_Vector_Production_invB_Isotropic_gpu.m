function dvec_y = FAME_Matrix_Vector_Production_invB_Isotropic_gpu(dvec_x, dB)
    dvec_y = dB.invB_eps.*dvec_x;
end