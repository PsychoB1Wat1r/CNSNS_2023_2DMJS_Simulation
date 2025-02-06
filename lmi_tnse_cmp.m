%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%% 数值求解代码 + 2023-9-4 %%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%% Yue-Yue Tao IEEE TNSE Proposition 1 Simulation %%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; % 清空命令行 工作区
tic  % 计时器
system_parameters_tnse % 导入系统参数
%% 循环定义决策变量
gamma = sdpvar(1, 1,'symmetric');   % 待求的性能指标gamma
% gamma = 0.4255^2;
Omega_h = sdpvar(1, 1,'symmetric');  Omega_v = sdpvar(1, 1,'symmetric');   
Omega = blkdiag(Omega_h, Omega_v);
Lh_g = sdpvar(1, 1,'symmetric');  Lv_g = sdpvar(1, 1,'symmetric');   
L = blkdiag(Lh_g, Lv_g);
for i = 1:NumSysMode
    for g = 1:NumConMode
        eval(['Ph_i',int2str(i),'= sdpvar(1, 1,''symmetric'');']);
        eval(['Pv_i',int2str(i),'= sdpvar(1, 1,''symmetric'');']);
        eval(['R_ig',int2str(i),int2str(g),'= sdpvar(2, 2,''symmetric'');']);
        eval(['K_g',int2str(g),'= sdpvar(1, 2,''full'');']);
    end
end
cnt = 0;  etch = cnt;  etcv = cnt;
beta = blkdiag(etch, etcv);  
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
        
        R_ig = eval(['R_ig',int2str(i), int2str(g)]);
        SumRig = eval(['R_ig_',int2str(i)]);
        K_g = eval(['K_g',int2str(g)]);
        %
        Phat = [sqrt(fiac(i,1))*P_i   sqrt(fiac(i,2))*P_i];
        
        hat_Fia_14 = [A*L + B*K_g  -B*K_g     E];      
        Pi_13 = [sqrt(fiao(i,1))*hat_Fia_14'   sqrt(fiao(i,2))*hat_Fia_14']';
        Pi_23 = [C*L + D*K_g  -D*K_g     F];
        Pi_33 = blkdiag(beta*Omega+R_ig-L-L', -Omega, -gamma*eye(1,1));
        % 线性矩阵不等式
        LMI_1 = [ -P_i         Phat;
                  (Phat)'   -SumRig];
        LMI_2 = [ -daigP        zeros(4,2)   Pi_13;
                   zeros(2,4)   -eye(2,2)    Pi_23;
                  (Pi_13)'      (Pi_23)'     Pi_33];
        Dim_LMI_1 = size(LMI_1);  Dim_LMI_2 = size(LMI_2);
        eval(['LMI_1_',int2str(i),int2str(g),'= [LMI_1<=delta*eye(Dim_LMI_1), LMI_2<=delta*eye(Dim_LMI_2), P_i>=0, R_ig>=0, Omega>=0];']); 
    end
end
%% 线性约束条件，需要同时满足
Constraint = [gamma>=0, LMI_1_11, LMI_1_12, LMI_1_21, LMI_1_22];
sdpsettings('solver', 'mosek','verbos',0); % 设置求解器为mosek，并打印少量信息
reuslt = optimize(Constraint, gamma); % 自动求解的优化函数
if reuslt.problem == 0 %problem = 0 代表求解成功
    [primal,~]=check(Constraint); % 检查约束
    if min(primal)>=0 && all(primal([1:2,4])>0) % 判断残差值
        disp('***********************************');
        disp('****Constraints are guaranteed*****');
        disp('**Primal residual is/are positive**');
        disp('***********************************');
        gamma = sqrt(double(gamma))
        K1 = value(K_g1)*inv(value(L)),  K2 = value(K_g2)*inv(value(L))
        Omega = value(Omega)
    else
        disp('********************************************');
        disp('**Warning: Primal residual is/are negative**');
        disp('********************************************');
        gamma = sqrt(double(gamma))
        K1 = value(K_g1)*inv(value(L)),  K2 = value(K_g2)*inv(value(L))
        Omega = value(Omega)
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