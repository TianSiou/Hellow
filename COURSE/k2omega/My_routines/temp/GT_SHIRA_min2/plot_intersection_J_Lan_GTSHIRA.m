
clear

theta = 0:0.001:2*pi;
x     = cos(theta)+1i*sin(theta);

load ../0804/EW_data_3p0
load near_c_data_eps_3p0 
load J_Lanczos_ew_no_intersection

idx = find(real(ew_unit_c(:,1)) <= 0.78 & real(ew_unit_c(:,1)) > 0.747);

figure(4) 
plot(ew,'bo','LineWidth',2);
hold on
plot(ew_unit_c(idx,1),'r*','LineWidth',1);
plot(tmp_ew,'g+','LineWidth',2);
plot(real(gamma),imag(gamma),'k^','LineWidth',1);
plot(x,'r-')
axis([0.45 0.78 0.62 0.9])
legend({'eig', 'GTSHIRA', 'J\_Lanczos', 'shift value'},'FontSize',12)
text(real(gamma)*1.005,imag(gamma)*1.005, ' \gamma','FontSize',14)
text(0.736,0.685, ' \eta_1','FontSize',14)
text(0.690,0.7264, ' \eta_2','FontSize',14)
text(0.661,0.7519, ' \eta_3','FontSize',14)
text(0.626,0.7818, ' \eta_4','FontSize',14)
text(0.582,0.815, ' \eta_5','FontSize',14)
text(0.557,0.7798, ' \eta_6','FontSize',14)
text(0.607,0.8492, ' \eta_7','FontSize',14)
text(0.459,0.8506, ' \eta_8','FontSize',14)
text(0.478,0.883, ' \eta_9','FontSize',14)

text(0.75,0.6638, ' \lambda_1','FontSize',14)
text(0.73,0.65, ' \lambda_2','FontSize',14)
text(0.75,0.63, ' \lambda_4','FontSize',14)

text(0.655,0.671, 'Lose to estimate \rightarrow','FontSize',12)
hold off

% open('ew_J_Lanc_GTSHIRA.fig')
% hold on
% 
% load J_Lanczos_ew_no_intersection
% plot(tmp_ew,'c+','LineWidth',2);
% plot(real(gamma),imag(gamma),'k^','LineWidth',2);
% axis([0.32 0.78 0.62 0.95])
% hold off
