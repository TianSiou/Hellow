function Test_Dense_Eigenproblems( mtx_gyroscopic, mtx_SHH, mtx_palindromic, mtx_symplecitc, mtx_ss, mtx_SH)
    EW_GQEP = polyeig( mtx_gyroscopic.K, mtx_gyroscopic.G, mtx_gyroscopic.M);
%     EW_SHH  = eig(full(mtx_SHH.H),full(mtx_SHH.S));
    EW_PQEP  = (1+EW_GQEP)./(1-EW_GQEP);
    EW_SS    = EW_PQEP + 1./EW_PQEP;
    
    subplot(1,3,1);
    plot(EW_GQEP,'b*');
    axis equal
    title('EW of GQEP(\mu)');
    subplot(1,3,2);
    plot(EW_PQEP,'b*');
    plot_unit_circle(1,3,2);
    axis equal
    title('EW of EW_PQEP matrix pair(\mu)');
    subplot(1,3,3);
    plot(EW_SS,'b*');
    plot_unit_circle(1,3,3);
    axis equal
    title('EW of EW_SS matrix pair(\eta)');
    
    
    
%     EW_PQEP = polyeig(mtx_palindromic.A, -mtx_palindromic.Q, mtx_palindromic.A.');
%     subplot(4,2,3);
%     plot(EW_PQEP,'b*');
%     plot_unit_circle(4,2,3);
%     axis equal
%     title('EW of PQEP(\nu)');
%     EW_GQEP_from_PQEP = -(1-EW_PQEP)./(1+EW_PQEP);
%     subplot(4,2,4);
%     plot(EW_GQEP_from_PQEP,'b*');
%     axis equal
%     title('EW of GQEP(\mu=-(1-\nu)/(1+\nu))');
%     
%     EW_Symp  = eig(full(mtx_symplecitc.M),full(mtx_symplecitc.L));
%     subplot(4,2,5);
%     plot(EW_Symp,'b*');
%     plot_unit_circle(4,2,5);
%     axis equal
%     title('EW of Symplectic(\nu)');
%     
%     EW_SS  = eig(full(mtx_ss.Ms),full(mtx_ss.Ls));
%     subplot(4,2,6);
%     plot(EW_SS,'b*');
%     plot_unit_circle(4,2,6);
%     axis equal
%     title('EW of Symplectic(\eta)');
%     EW_PQEP_from_SS = 0.5*(EW_SS-sqrt(EW_SS.^2-4));
%     EW_PQEP_from_SS = [EW_PQEP_from_SS;1./EW_PQEP_from_SS];
%     subplot(4,2,7);
%     plot(EW_PQEP_from_SS,'b*');
%     plot_unit_circle(4,2,7);
%     axis equal
%     title('EW of PQEP(\nu = 0.5*(\eta-(\eta^2-4)^{0.5}) and 0.5*(\eta+(\eta^2-4)^{0.5}) )');
%     
%     EW_SHpair  = eig(full(mtx_SH.K),full(mtx_SH.N));
%     subplot(4,2,8);
%     plot(EW_SHpair,'b*');
%     plot_unit_circle(4,2,8);
%     axis equal
%     title('EW of skew-Hamiltonian pair(\eta)');
    
%     nu_0 = -0.5;
%     eta_0 = nu_0 + 1/nu_0;
%     EW_SHpair_shift  = eig(full(mtx_SH.K_hat(nu_0)),full(mtx_SH.N_hat(nu_0)));
%     EW_SHpair_shift  = 1./EW_SHpair_shift + eta_0;
%     subplot(4,2,8); hold on;
%     plot(EW_SHpair_shift,'ro');
%     title('EW of skew-Hamiltonian pair(\eta)');
pause(0.1);
end
function plot_unit_circle(a,b,c)
theta = 0:2*pi/100:2*pi;
x = cos(theta);
y = sin(theta);

subplot(a,b,c); hold on;
plot(x,y,'k-');

    
end
    