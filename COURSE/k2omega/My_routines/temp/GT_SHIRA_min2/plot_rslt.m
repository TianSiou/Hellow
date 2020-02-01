clear
% load all_eps_unit_ew.mat
load all_eps_unit_ew_defl_20 %all_eps_unit_ew_0715_all.mat

figure(1)
plot(Eps_interval,no_unit_ew,'bo-','LineWidth',2);
ylabel('Numbers of unimodular eigenvalues','FontSize',14)
xlabel('$\mathcal{E}$','Interpreter','LaTex', 'FontSize', 14)

% h=legend({'$a_i$', '$b_i$'},'FontSize',14,'Location','East');
% set(h,'Interpreter','latex')

figure(2)
plot(Eps_interval,iter_all,'bo-','LineWidth',2);
xlabel('$\mathcal{E}$','Interpreter','LaTex', 'FontSize', 14)
ylabel('Numbers of using GTSHIRA','FontSize',14)
%set(h,'Interpreter','latex')

figure(3)
plot(Eps_interval,ratio_deuplicate_ew,'bo','LineWidth',2);
xlabel('$\mathcal{E}$','Interpreter','LaTex', 'FontSize', 14)
ylabel('Percentage of deuplicated eigenvalues','FontSize',14)

figure(4)
plot(Eps_interval,CPUTime_all,'bo','LineWidth',2);
xlabel('$\mathcal{E}$','Interpreter','LaTex', 'FontSize', 14)
ylabel('CPU time $T_d$ (sec.)','Interpreter','LaTex','FontSize',14)

CPUTime_all_defl = CPUTime_all;

clear CPUTime_all Eps_interval iter_all ratio_deuplicate_ew

load all_eps_unit_ew_without_defl

diff_CPUTime = CPUTime_all_defl - CPUTime_all;
figure(5)
plot(Eps_interval,diff_CPUTime,'bo','LineWidth',2);
xlabel('$\mathcal{E}$','Interpreter','LaTex', 'FontSize', 14)
ylabel('$T_d - T$','Interpreter','LaTex','FontSize',14)

figure(2)
hold on
plot(Eps_interval,iter_all,'rx-','LineWidth',2);
legend({'with deflation', 'without deflation'},'FontSize',14,'Location','NorthWest');
hold off

figure(3)
hold on 
plot(Eps_interval,ratio_deuplicate_ew,'rx','LineWidth',2);
legend({'with deflation', 'without deflation'},'FontSize',14,'Location','SouthWest');
hold off

figure(4)
% hold on
% plot(Eps_interval,CPUTime_all,'rx','LineWidth',2);
hold off