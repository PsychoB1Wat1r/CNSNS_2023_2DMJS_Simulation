%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%% 双模态三维性能指标 gamma 绘制 + 2023-9-4 %%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; % 清空命令行 工作区
tic  % 计时器
system_parameters_tnse % 导入系统参数
%% 三维性能指标绘图参数
T = 1; dt = 0.25; Nt = T/dt; % 坐标轴1
Lx = 1; dx = 0.25; Nx = Lx/dx; % 坐标轴2
%% （异步） 马尔可夫过程alpha 系统切换 跳变参数参数一 
Phi = [0.9 0.1;  0.8 0.2];
%% 跳变参数二
beta11 = 0:dt:T; beta12 = 1-beta11; % 切分坐标轴
beta21 = 0:dx:Lx; beta22 = 1-beta21; % 切分坐标轴
%% 丢包率
k = 0.8; hatk = sqrt(k*(1-k));
%% 决策变量
gamma = sdpvar(1, 1,'full'); % 待求的性能指标gamma
%% 区间切分
fia11 = zeros(Nt+1,Nx+1); fia12 = zeros(Nt+1,Nx+1);
fia21 = zeros(Nt+1,Nx+1); fia22 = zeros(Nt+1,Nx+1);
gammaset = zeros(Nt+1,Nx+1); % 用于存储所计算出的性能指标Gamma值的矩阵 便于绘图
for i = 1:Nt+1
    for j = 1:Nx+1
    fia11(i,j) = beta11(j); fia12(i,j) = beta12(j);
    fia21(i,j) = beta21(i); fia22(i,j) = beta22(i);
    end
end
%% 系统模态/控制器模态
NumSysMode = 2; NumConMode = 2;
for i = 1:NumSysMode
    for g = 1:NumConMode
        eval(['Ph_i',int2str(i),'= sdpvar(1, 1,''symmetric'');']);
        eval(['Pv_i',int2str(i),'= sdpvar(1, 1,''symmetric'');']);
        eval(['R_ig',int2str(i),int2str(g),'= sdpvar(2, 2,''symmetric'');']);
        eval(['Lh_g',int2str(g),'= sdpvar(1, 1,''symmetric'');']);
        eval(['Lv_g',int2str(g),'= sdpvar(1, 1,''symmetric'');']);
        eval(['K_g',int2str(g),'= sdpvar(1, 2,''full'');']);
        eval(['Omega_h',int2str(i),'= sdpvar(1, 1,''symmetric'');']);
        eval(['Omega_v',int2str(i),'= sdpvar(1, 1,''symmetric'');']);
    end
end
etc = 0.3; % 事件触发参数
beta = blkdiag(etc, etc);  
%
R_ig_1 = blkdiag(R_ig11, R_ig12);  R_ig_2 = blkdiag(R_ig21, R_ig22);
%
P_i_1 = blkdiag(Ph_i1, Pv_i1);  P_i_2 = blkdiag(Ph_i2, Pv_i2);
daigP = blkdiag(P_i_1, P_i_2);
%% 循环定义LMI
for i = 1:NumSysMode
    for g = 1:NumConMode
        A = eval(['A',int2str(i)]);  B = eval(['B',int2str(i)]);  E = eval(['E',int2str(i)]);
        C = eval(['C',int2str(i)]);  D = eval(['D',int2str(i)]);  F = eval(['F',int2str(i)]);
        Ph_i = eval(['Ph_i',int2str(i)]);  Pv_i = eval(['Pv_i',int2str(i)]);
        P_i = blkdiag(Ph_i, Pv_i);
        
        Omega_h = eval(['Omega_h',int2str(i)]);  Omega_v = eval(['Omega_v',int2str(i)]);
        Omega = blkdiag(Omega_h, Omega_v);
        
        Lh_g = eval(['Lh_g',int2str(g)]);  Lv_g = eval(['Lv_g',int2str(g)]);
        L_g = blkdiag(Lh_g, Lv_g);
        
        R_ig = eval(['R_ig',int2str(i), int2str(g)]);
        SumRig = eval(['R_ig_',int2str(i)]);
        K_g = eval(['K_g',int2str(g)]);

        hat_Fia_14 = [A*L_g + k*B*K_g  -k*B*K_g     E];
        hat_Fia_24 = [hatk*B*K_g       -hatk*B*K_g  zeros(2,1)];         
        Pi_14 = [sqrt(Phi(i,1))*hat_Fia_14'   sqrt(Phi(i,2))*hat_Fia_14']';
        Pi_24 = [sqrt(Phi(i,1))*hat_Fia_24'   sqrt(Phi(i,2))*hat_Fia_24']';
        Pi_34 = [C*L_g + k*D*K_g  -k*D*K_g     F;
                 hatk*D*K_g       -hatk*D*K_g  zeros(2,1)];
        Pi_44 = blkdiag(beta*Omega+R_ig-L_g'-L_g', -Omega, -gamma*eye(1,1));
        LMI_2 = [ -daigP        zeros(4,4)   zeros(4,4)  Pi_14;
                   zeros(4,4)  -daigP        zeros(4,4)  Pi_24;
                   zeros(4,4)   zeros(4,4)  -eye(4,4)    Pi_34;
                  (Pi_14)'     (Pi_24)'     (Pi_34)'     Pi_44];
        eval(['Condition',int2str(i),int2str(g),'= [LMI_2<=0, P_i>=0, R_ig>=0, Omega>=0];']); 
    end
end
%% 重复迭代计算gamma的过程（二重循环迭代）
for i = 1:Nt+1
    for j = 1:Nx+1
        Tp_1 = [sqrt(fia11(i,j))*P_i_1   sqrt(fia12(i,j))*P_i_1];
        Tp_2 = [sqrt(fia21(i,j))*P_i_2   sqrt(fia22(i,j))*P_i_2];
        LMI1_1 = [-P_i_1     Tp_1;
                   Tp_1'    -R_ig_1];
        LMI1_2 = [-P_i_2     Tp_2;
                   Tp_2'    -R_ig_2];
       %% 约束条件求解，需要同时满足
        Constraints = [gamma>=0, LMI1_1<=0, LMI1_2<=0, Condition11,  Condition12, Condition21,  Condition22];
        sdpsettings('solver', 'mosek','verbos',0); % 设置求解器为mosek
        reuslt = optimize(Constraints, gamma); % 自动求解的函数
        if reuslt.problem == 0 % problem = 0 代表求解成功
            gammaset(i,j) = sqrt(double(gamma));  % 将计算出的Gamma值存储在矩阵gammaset中，方便绘图
        else
        end
    end
end
%% 性能指标gamma 三维图
H1 = figure;
[x,y] = meshgrid(0:dt:T,0:dx:Lx);
z = gammaset(1:Nt+1,1:Nx+1);
mesh(x, y, z, 'LineWidth', 1, 'marker', '.', 'MarkerFaceColor',[0.5,0.5,0.5], 'markersize', 10);
xlabel('$\pi_{21}$','Interpreter','latex','Fontsize',19);
ylabel('$\pi_{11}$','Interpreter','latex','Fontsize',19);
zlabel('$\gamma_{c}^{*}$','Interpreter','latex','Fontsize',19);
axis([0 T 0 Lx 1.52 1.83]);  grid on;  grid minor;  box on;
colormap(hsv(3));
%% 性能指标gamma 二维图
H2 = figure;
plot((0:dt:T), diag(gammaset(1:Nt+1,1:Nx+1)), 'color', '#083AA9', 'linestyle', '-.', 'linewidth', 2);
hold on
plot((0:dt:T), diag(rot90(gammaset(1:Nt+1,1:Nx+1))), 'color', '#FF2E63', 'linestyle', '-.', 'linewidth', 2, 'marker', 'p', 'MarkerFaceColor',[0.5,0.5,0.5], 'markersize', 4.5);
hold on
plot([0 T],[1.53956 1.53956], 'color', '#222831', 'linewidth', 2);
set(gca,'FontSize',11, 'Linewidth', 1);
xlabel('$\pi_{11}$','Interpreter','latex','Fontsize', 19);
ylabel('$\gamma_{c}^{*}$','Interpreter','latex','Fontsize', 19);
h1 = legend('Mode-independent $(\pi_{11}=\pi_{21})$', 'Asynchronous $(\pi_{11}+\pi_{21}=1)$', 'Synchronous $(\pi_{11}=1$ and $\pi_{21}=0)$', 'Interpreter','latex','Fontsize',13, 'location', 'best');
axis([0 T 1.52 1.82]);  grid on;  grid minor;  box on;
disp(['运行时间: ',num2str(toc)]); 