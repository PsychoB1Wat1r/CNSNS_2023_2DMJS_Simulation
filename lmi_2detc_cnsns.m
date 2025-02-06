%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%% 数值求解代码 + 2023-9-4 %%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%% 2D MJS + 事件触发 + 丢包 + 异步控制 %%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; % 清空命令行 工作区
tic  % 计时器
system_parameters_tnse % 导入系统参数
%% 循环定义决策变量
% gamma = sdpvar(1, 1,'symmetric');   % 待求的性能指标gamma
gamma = 1.7834^2; % 给定gamma值 用于求解控制增益 与 触发矩阵
for i = 1:NumSysMode
    for g = 1:NumConMode
        eval(['Ph_i',int2str(i),'= sdpvar(1, 1,''symmetric'');']);
        eval(['Pv_i',int2str(i),'= sdpvar(1, 1,''symmetric'');']);
        eval(['R_ig',int2str(i),int2str(g),'= sdpvar(2, 2,''symmetric'');']);
        eval(['Lh_g',int2str(g),'= sdpvar(1, 1,''symmetric'');']);
        eval(['Lv_g',int2str(g),'= sdpvar(1, 1,''symmetric'');']);
        eval(['K_g',int2str(g),'= sdpvar(1, 2,''full'');']);
        eval(['Omega_h',int2str(i),'= sdpvar(1, 1,''symmetric'');']); % 模态相关事件触发权值矩阵
        eval(['Omega_v',int2str(i),'= sdpvar(1, 1,''symmetric'');']);
    end
end
% 中间变量
R_ig_1 = blkdiag(R_ig11, R_ig12);  R_ig_2 = blkdiag(R_ig21, R_ig22);
% 中间变量
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
        %
        Phat = [sqrt(fiac(i,1))*P_i   sqrt(fiac(i,2))*P_i];
        
        hat_Fia_14 = [A*L_g + k*B*K_g  -k*B*K_g     E];
        hat_Fia_24 = [hatk*B*K_g       -hatk*B*K_g  zeros(2,1)];         
        Pi_14 = [sqrt(fiao(i,1))*hat_Fia_14'   sqrt(fiao(i,2))*hat_Fia_14']';
        Pi_24 = [sqrt(fiao(i,1))*hat_Fia_24'   sqrt(fiao(i,2))*hat_Fia_24']';
        Pi_34 = [C*L_g + k*D*K_g  -k*D*K_g     F;
                 hatk*D*K_g       -hatk*D*K_g  zeros(2,1)];
        Pi_44 = blkdiag(etc*Omega+R_ig-L_g-L_g', -Omega, -gamma*eye(1,1));
        % 线性矩阵不等式
        LMI_1 = [ -P_i         Phat;
                  (Phat)'   -SumRig];
        LMI_2 = [ -daigP        zeros(4,4)   zeros(4,4)  Pi_14;
                   zeros(4,4)  -daigP        zeros(4,4)  Pi_24;
                   zeros(4,4)   zeros(4,4)  -eye(4,4)    Pi_34;
                  (Pi_14)'     (Pi_24)'     (Pi_34)'     Pi_44];
        eval(['LMI_1_',int2str(i),int2str(g),'= [LMI_1<=delta*eye(6), LMI_2<=delta*eye(17), P_i>=0, R_ig>=0, Omega>=0];']); 
    end
end
%% 线性约束条件，需要同时满足
Constraint = [gamma>=0, LMI_1_11, LMI_1_12, LMI_1_21, LMI_1_22];
sdpsettings('solver', 'mosek', 'verbos', 0); % 设置求解器为mosek，并打印少量信息
reuslt = optimize(Constraint, gamma); % 自动求解的优化函数
if reuslt.problem == 0 %problem = 0 代表求解成功
    [primal,~]=check(Constraint); % 检查约束
    if min(primal)>=0 && all(primal([1:2,4])>0) % 判断残差值
        disp('***********************************');
        disp('****Constraints are guaranteed*****');
        disp('**Primal residual is/are positive**');
        disp('***********************************');
        gamma = sqrt(double(gamma))
        K1 = value(K_g1)*inv(value(blkdiag(Lh_g1, Lv_g1))),  K2 = value(K_g2)*inv(value(blkdiag(Lh_g2, Lv_g2)))
        Omega1 = inv(value(blkdiag(Lh_g1, Lv_g1)'))*value(blkdiag(Omega_h1, Omega_v1))*inv(value(blkdiag(Lh_g1, Lv_g1)))
        Omega2 = inv(value(blkdiag(Lh_g2, Lv_g2)'))*value(blkdiag(Omega_h2, Omega_v2))*inv(value(blkdiag(Lh_g2, Lv_g2)))
    else
        disp('********************************************');
        disp('**Warning: Primal residual is/are negative**');
        disp('********************************************');
        gamma = sqrt(double(gamma))
        K1 = value(K_g1)*inv(value(blkdiag(Lh_g1, Lv_g1))),  K2 = value(K_g2)*inv(value(blkdiag(Lh_g2, Lv_g2)))
        Omega1 = inv(value(blkdiag(Lh_g1, Lv_g1)'))*value(blkdiag(Omega_h1, Omega_v1))*inv(value(blkdiag(Lh_g1, Lv_g1)))
        Omega2 = inv(value(blkdiag(Lh_g2, Lv_g2)'))*value(blkdiag(Omega_h2, Omega_v2))*inv(value(blkdiag(Lh_g2, Lv_g2)))   
        check(Constraint); %检查残余值
    end
else
    disp('**************************************');
    disp('****Constraints are not guaranteed****');
    disp('**************************************'); 
    reuslt.info
    yalmiperror(reuslt.problem) % 打印出错信息
end
disp(['运行时间: ',num2str(toc)]); 