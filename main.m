clc;
clear;
y0=[0 0 1 0 0 0 0 0 0];
tspan=[0,12];
tau=[1 2 3 2];%time lag of y(2) y(5) y(7) y(9)
sol=dde23(@test1,tau,y0,tspan);
T = linspace(tspan(1),tspan(2),1000);
Y = deval(sol,T);
figure;
plot(T,Y,'LineWidth',1);
    xlabel('Time /s');
    ylabel('miumol/L');
    legend('mRNA_H_L_R','pt_H_L_R','DNA_C_1_2_0_S','mRNA_d_e_l_t_a','pt_d_e_l_t_a','mRNA_R_e_l_B','pt_R_e_l_B','mRNA_R_e_l_E','pt_R_e_l_E','Location','NorthWest');
    hold on;
grid on
%set(gca,'xtick',0:1:tspan(1,2),'ytick',0:.5:3.5)
hold off


