clc;
clear;
%(1)mRNA_HLR
q_mRNA_HLR=0.1;
C_mRNA_HLR=0.3;
%(2)pt_HLR
r_pt_HLR=0.2;
q_HLR_pt=0.01;
k_deC=0.1;

for Ls=1:100;

%timeLapse todo 
%(3)DNA_C120S
%(4)mRNA_delta
r_mRNA_delta_1=0.2;
r_mRNA_delta_2=0.25;
q_mRNA_delta=0.01;
%(5)pt_delta
r_pt_delta=0.2;
q_pt_delta=0.1;
%(6)mRNA_RelE
r_mRNA_delta_2=0.2;
q_mRNA_RelB=0.1;
%(7)pt_RelB
r_pt_RelB=0.2;
q_pt_RelB=0.01;
kd=0.1;
%(8)mRNA_RelE
C_mRNA_RelE=0.3;
q_mRNA_RelE=0.1;
%(9)pt_RelE
r_mRNA_RelE=0.1;
q_pt_RelE=0.1;

sig = [0;0.1;0;0;0;0;0;0;0];% ??
dt = 1e-2;
t = 0:dt:10;
dtt = 0:1e-3:1;
f = @(t,y)[C_mRNA_HLR-q_mRNA_HLR*y(1);%1    
    r_pt_HLR*y(1)-q_HLR_pt*y(2)-k_deC*y(3)*y(2)^2+L_repress(Ls)*(1-y(3));%2
    -k_deC*y(3)*y(2)^2+L_repress(Ls)*(1-y(3));%3
    r_mRNA_delta_1*y(3)+r_mRNA_delta_2*y(5)-q_mRNA_delta*y(4);%4
    y(4)*r_pt_delta-q_pt_delta*y(5);%5  
    r_mRNA_delta_2*y(5)-q_mRNA_RelB*y(6);%6
    r_pt_RelB*y(6)-q_pt_RelB*y(7)-kd*y(7)*y(9);%7
    C_mRNA_RelE-q_mRNA_RelE*y(8);%8
    y(8)*r_mRNA_RelE-y(9)*q_pt_RelE-kd*y(9)*y(7)%9
    ];
g = @(t,y)sig.*y;
opts = sdeset('DGFUN',sig);
y = sde_milstein(f,g,t,[0 0 1 0 0 0 0 0 0]);
%figure;
% for i=1:9
%     if i==7
%         plot(t,y(:,i),'LineWidth',2);
%     else
%         plot(t,y(:,i),'LineWidth',1);
%     end
%     xlabel('Time /s');
%     ylabel('miumol/L');
%     legend('mRNA_H_L_R','pt_H_L_R','DNA_C_1_2_0_S','mRNA_d_e_l_t_a','pt_d_e_l_t_a','mRNA_R_e_l_B','pt_R_e_l_B','mRNA_R_e_l_E','pt_R_e_l_E','Location','NorthWest');
%     hold on;
% end
%plot(t,y);
a(1,Ls)=y(800,7);
end