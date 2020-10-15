function dy=test1(t,y,Z)
%time lag
ylag1 = Z(:,1);
ylag2 = Z(:,2);
ylag3 = Z(:,3);
ylag4 = Z(:,4);
%(1)mRNA_HLR
q_mRNA_HLR=0.1;
C_mRNA_HLR=0.3;
%(2)pt_HLR
r_pt_HLR=0.2;
q_HLR_pt=0.01;
k_deC=0.1;
Ls=1;
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
dy=[C_mRNA_HLR-q_mRNA_HLR*y(1);
    r_pt_HLR*y(1)-q_HLR_pt*ylag1(2)-k_deC*y(3)*ylag1(2)^2+L_repress(Ls)*(1-y(3))*randn(1);
    -k_deC*y(3)*ylag1(2)^2+L_repress(Ls)*(1-y(3))*randn(1);
    r_mRNA_delta_1*y(3)+r_mRNA_delta_2*ylag2(5)-q_mRNA_delta*y(4);
    y(4)*r_pt_delta-q_pt_delta*ylag2(5);
    r_mRNA_delta_2*ylag2(5)-q_mRNA_RelB*y(6);
    r_pt_RelB*y(6)-q_pt_RelB*ylag3(7)-kd*ylag3(7)*ylag4(9);
    C_mRNA_RelE-q_mRNA_RelE*y(8);
    y(8)*r_mRNA_RelE-ylag4(9)*q_pt_RelE-kd*ylag4(9)*ylag3(7);
    ];
end
% translation speed = almost 15 aa/s * 
% binding ribosome number * 
% available mRNA number
