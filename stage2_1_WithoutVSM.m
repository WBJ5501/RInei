clear all
clc

run ieee_33_node_system.m;
run DG_Load_MPC.m;
run oltc_cbs_and_withoutvsm.m;
load('风光负荷_96_112_15min_随机性.mat');

tic
t1 = clock;
%% 常量定义
%------------------ 基准值设置--------------------------------------
U_b=12.66;%kV
S_b=1;%MW
Z_b=U_b^2/S_b;
%----------------电阻电抗--------------------------------
r_ij = Branch( : , 4 )/Z_b;
x_ij = Branch( : , 5 )/Z_b;
%------------------目标函数值--------------------------------------
N =8;%时间
result = zeros(1,N);
%------------------电压上下限--------------------------------------
a1=0.95;
a2=1.05;
%----------------变压器抽头及电容器--------------------------------
delta_T = 1;%hour
K_max = 5;%抽头最大挡位
kk = (0:K_max*2);
k_ij0 = 1;%初始变比

delta_kij = 0.1/(K_max*2);%每个挡位电压变化率
%------------------电容器单位无功功率--------------------------------------
qCB = 0.06;%每个电容器的无功功率（Mvar）

%----------------光伏安装位置及其容量（Mwh）--------------------------------
q_n=6;

pv1=10;
pv2=14;
pv3=17;
pv4=24;
pv5=28;%光伏安装位置
pv6=32;

S_solar=zeros(1,q_n);
S_solar(1)=0.8*1.1;%光伏容量
S_solar(2)=0.3*1.1;
S_solar(3)=0.3*1.1;
S_solar(4)=2.3*1.1;
S_solar(5)=3.5*1.1;
S_solar(6)=1*1.1;
%----------------光伏无功出力限制值--------------------------------
Q_pv=zeros(1,q_n);
for q_pv=1:q_n
    Q_pv(q_pv)=S_solar(q_pv)*0.436*0.9;
end
%----------------SOP无功出力限制值--------------------------------
% Q_sop = 2*0.6*0.9;
%---------------------ESSOP各子系统容量和损耗系数-------------------------------
%------容量------------
S_acdc = 2;
S_dcdc = 0.6;
S_ess = 1;
%------损耗系数------------
C_a0=0.005;
C_a1=0.005;
C_a2=0.02;

C_d0=0.0035;
C_d1=0.0025;
C_d2=0.014;

% C_aux=0.002*S_ess;
C_aux=0;
C_ES=0.016;

S_ESS = zeros(96,N+1);
%% 多时间尺度循环
for t=1:96

 Ktij = zeros(1,N);
 NtCB =zeros(3,N);
  
for c = t:(t+N-1)
    Ktij( 1 , c-t+1 ) = Ktij_1( c ); 
    
    NtCB( 1 , c-t+1 ) = NtCB_1( c ); 
    NtCB( 2 , c-t+1 ) = NtCB_2( c );
    NtCB( 3 , c-t+1 ) = NtCB_3( c ); 
end 

%----------------负荷光伏数据及恒定功率因数--------------------------------
p_Solar = zeros(33,N);
p_Wind  = zeros(33,N);
p_Load  = zeros(33,N);
q_Load  = zeros(33,N);
tan_theta = q_load ./ p_load;%恒定功率因数
for b = t:(t+N-1)
    p_Solar( : , b-t+1 ) = Solar_radio * Solar_data(t, b );
    p_Load ( : , b-t+1 ) = Load_radio  *  Load_data(t, b );
    q_Load ( : , b-t+1 ) = ( p_Load ( : , b-t+1 ) ./ p_load ) .* ( q_load );
end

%% 定义初始运行点决策变量
q_Solar =sdpvar(q_n, N, 'full');
%-----------------光伏削减量------------------------------------
PV_cut=sdpvar(q_n, N, 'full');
P_pv=sdpvar(q_n, N, 'full');
%----------------电压、支路有功、支路无功--------------------------------
x_ui_square = sdpvar(33, N, 'full');
x_pij = sdpvar(32, N, 'full');
x_qij = sdpvar(32, N, 'full');
%----------------变压器抽头--------------------------------
% btij = binvar(K_max*2+1, N, 'full');
% 
% Ktij = intvar(1, N+1, 'full');
% Ktij1 = intvar(1, N, 'full');
% Ktij2 = intvar(1, N, 'full');
%----------------电容器--------------------------------
% qc=6;
% NtCB = intvar(qc, N+1, 'full');
% NtCB1 = intvar(qc, N, 'full');
% NtCB2 = intvar(qc, N, 'full');
%--------------------1组ESSOP/双馈线-------------------------------------
sv=2;
%VSC接通/断开变量
b1 = binvar(sv, N, 'full');
b2 = binvar(sv, N, 'full');
b3 = binvar(sv-1, N, 'full');
%AC/DC视在功率辅助变量/非负
k11= sdpvar(sv, N, 'full');
k12= sdpvar(sv, N, 'full');
k21= sdpvar(sv, N, 'full');
k22= sdpvar(sv, N, 'full');
%四个子系统损耗
P1_L=sdpvar(sv, N, 'full');
P2_L=sdpvar(sv, N, 'full');
P3_L=sdpvar(sv-1, N, 'full');
P4_L=sdpvar(sv-1, N, 'full');
%四个子系统有功功率
P1=sdpvar(sv, N, 'full');
P2=sdpvar(sv, N, 'full');
P3=sdpvar(sv-1, N, 'full');
P4=sdpvar(sv-1, N, 'full');
%AC/DC无功功率
Q1=sdpvar(sv, N, 'full');
Q2=sdpvar(sv, N, 'full');
%DC/DC损耗中的绝对值项和平方项
P3_abs=sdpvar(sv-1, N, 'full');
k3=sdpvar(sv-1, N, 'full');
%ESS损耗中的平方项和荷电状态
k4=sdpvar(sv-1, N, 'full');    
Soc_ess=sdpvar(sv-1, N+1, 'full');
%-----------------电压偏差------------------------------------
Aux = sdpvar(33, N, 'full');
%----------------损耗-----------------------------
P_loss=sdpvar(32, N, 'full'); 

%% 定义约束条件
Constraints = [];
% Constraints = [ Constraints , q_Solar== 0 ];
% Constraints = [ Constraints , b1 == 0 ];
% Constraints = [ Constraints , b2 == 0 ];
% Constraints = [ Constraints , b3 == 0 ];
for opt_num = 1: (N)
     
    tic
%% 目标函数


        f(opt_num)=1000*0.25*(0.08*( sum(P_loss( : , opt_num ))+sum( P1_L( : , opt_num )+ P2_L( : , opt_num ) )+ P3_L( : , opt_num )+ P4_L( : , opt_num ))...
                  + 0.64*sum(PV_cut( : , opt_num )))+25*(sum( Aux( : , opt_num ))); 
              
        f_start(opt_num)= sum(P_loss( : , opt_num ))*0.25;
        f_essop(opt_num)= (sum( P1_L( : , opt_num )+ P2_L( : , opt_num ))+ P3_L( : , opt_num )+ P4_L( : , opt_num ))*0.25; 
        f_PV_cut(opt_num)= sum(PV_cut( : , opt_num ))*0.25;

        
%% 平衡节点电压约束
   Constraints = [ Constraints , x_ui_square( 1 , opt_num ) == 1 ];
   
%% 有功无功平衡
for k = 2 : 33
    node_out = find(Branch(:,2) == k);%下游支路
    node_in  = find(Branch(:,3) == k);%上游支路

        Constraints = [ Constraints ,P_loss(k-1 , opt_num ) >= r_ij( k-1 )*(x_pij( k-1 , opt_num )^2 + x_qij( k-1 , opt_num )^2) ]; 
    
    if( k == 2 )
        Constraints = [ Constraints ,x_pij( node_in , opt_num ) + p_Solar( k , opt_num ) + p_Wind( k , opt_num )...
                                     - p_Load( k , opt_num ) == sum( x_pij( node_out , opt_num ) ) ];
        Constraints = [ Constraints , x_qij( 1 , opt_num )...
                                      - q_Load( k, opt_num  )== sum( x_qij( node_out , opt_num ) ) ];
                                  
                           
%% ESSOP1    
    elseif( k == 12 )
        Constraints = [ Constraints , x_pij( node_in , opt_num )  ...
                                      + p_Solar( k , opt_num ) + p_Wind( k , opt_num )...
                                      + P1( 1 , opt_num )- p_Load( k , opt_num ) == sum( x_pij( node_out , opt_num ) ) ];
        Constraints = [ Constraints , x_qij( node_in ,opt_num )  ...
                                      + Q1( 1 , opt_num )- q_Load( k , opt_num  ) == sum( x_qij( node_out , opt_num ) ) ];
      
    
    elseif( k == 22 )
        Constraints = [ Constraints , sum( x_pij( node_in ,opt_num ) )...
                                       + P2( 1 , opt_num )+ p_Solar( k , opt_num ) + p_Wind( k , opt_num ) - p_Load( k , opt_num ) == 0 ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num )  )...
                                       + Q2( 1 , opt_num )- q_Load( k , opt_num  )== 0 ];
                     
%% ESSOP2 
    elseif( k == 25 )
        Constraints = [ Constraints , sum( x_pij( node_in , opt_num )  )...
                                     + P1(2 , opt_num )+ p_Solar( k , opt_num ) + p_Wind( k , opt_num ) - p_Load( k , opt_num ) == 0 ];
        Constraints = [ Constraints , sum( x_qij( node_in ,opt_num )  )...
                                     + Q1( 2 , opt_num )- q_Load( k , opt_num ) == 0 ];                  
                                
    elseif( k == 29 )
        Constraints = [ Constraints , sum( x_pij( node_in , opt_num ) )...
                                     + P2( 2 , opt_num )+ p_Solar( k , opt_num ) + p_Wind( k , opt_num )...
                                    - p_Load( k , opt_num ) == sum( x_pij( node_out , opt_num ) ) ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num ) )...
                                     + Q2( 2 , opt_num )- q_Load( k , opt_num )== sum( x_qij( node_out , opt_num ) ) ];                     
                                
    elseif( k == 18 )
        Constraints = [ Constraints , sum( x_pij( node_in , opt_num )  )...
                                    + p_Solar( k , opt_num ) + p_Wind( k , opt_num ) - p_Load( k , opt_num ) == 0 ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num )  )...
                                    - q_Load( k , opt_num ) == 0 ];                        
  
    elseif( k == 33 )
        Constraints = [ Constraints , sum( x_pij( node_in , opt_num ) )...
                                    + p_Solar( k , opt_num ) + p_Wind( k , opt_num ) - p_Load( k , opt_num ) == 0 ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num ) )...
                                    - q_Load( k , opt_num  ) == 0 ];

                                                                                                                
    elseif( k == 6 )
        Constraints = [ Constraints , sum( x_pij( node_in , opt_num )  )...
                                      + p_Wind( k , opt_num ) - p_Load( k , opt_num ) == sum( x_pij( node_out , opt_num ) ) ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num )  )...
                                      - q_Load(k , opt_num )+ NtCB(1,opt_num  )*qCB   == sum( x_qij( node_out , opt_num ) ) ];
                                  
    elseif( k == 13 )
        Constraints = [ Constraints , sum( x_pij( node_in , opt_num )  )...
                                      + p_Solar( k , opt_num ) + p_Wind( k , opt_num ) - p_Load( k , opt_num ) == sum( x_pij( node_out , opt_num ) ) ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num ) )...
                                      - q_Load(k , opt_num )+ NtCB( 2,opt_num   )*qCB   == sum( x_qij( node_out , opt_num ) ) ];
                           

    elseif( k == 31 )
        Constraints = [ Constraints , sum( x_pij( node_in , opt_num )  )...
                                      + p_Solar( k , opt_num ) + p_Wind( k , opt_num ) - p_Load( k , opt_num ) == sum( x_pij( node_out , opt_num ) ) ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num )  )...
                                      - q_Load(k , opt_num )  + NtCB(3, opt_num   )*qCB  == sum( x_qij( node_out , opt_num ) ) ];
                                                      
                                         
%% 含风光机组的节点有功无功功率平衡    
    elseif( k == pv1 )
        Constraints = [ Constraints , sum( x_pij( node_in , opt_num ) )...
                                      + p_Solar( k , opt_num )-PV_cut( 1 , opt_num ) + p_Wind( k , opt_num ) - p_Load( k , opt_num ) == sum( x_pij( node_out , opt_num ) ) ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num )  )...
                                      - q_Load(k , opt_num )+ q_Solar(1 , opt_num ) == sum( x_qij( node_out , opt_num ) ) ];
                         
                                  
    elseif( k ==pv2 )
        Constraints = [ Constraints , sum( x_pij( node_in , opt_num )  )...
                                      + p_Solar( k , opt_num )-PV_cut( 2 , opt_num ) + p_Wind( k , opt_num ) - p_Load( k , opt_num ) == sum( x_pij( node_out , opt_num ) ) ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num )  )...
                                      - q_Load(k , opt_num )+ q_Solar(2 , opt_num ) == sum( x_qij( node_out , opt_num ) ) ];

    elseif( k ==pv3 )
        Constraints = [ Constraints , sum( x_pij( node_in , opt_num )  )...
                                      + p_Solar( k , opt_num )-PV_cut( 3 , opt_num ) + p_Wind( k , opt_num ) - p_Load( k , opt_num ) == sum( x_pij( node_out , opt_num ) ) ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num )  )...
                                      - q_Load(k , opt_num )+ q_Solar( 3 , opt_num ) == sum( x_qij( node_out , opt_num ) ) ];    
                                                         
    elseif( k ==pv4)
        Constraints = [ Constraints , sum( x_pij( node_in , opt_num )  )...
                                      + p_Solar( k , opt_num )-PV_cut( 4 , opt_num ) + p_Wind( k , opt_num ) - p_Load( k , opt_num ) == sum( x_pij( node_out , opt_num ) ) ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num ) )...
                                      - q_Load(k , opt_num )+ q_Solar(4 , opt_num ) == sum( x_qij( node_out , opt_num ) ) ];
                                                                                            
     elseif( k == pv5 )
        Constraints = [ Constraints , sum( x_pij( node_in , opt_num ) )...
                                      + p_Solar( k , opt_num )-PV_cut( 5 , opt_num ) + p_Wind( k , opt_num ) - p_Load( k , opt_num ) == sum( x_pij( node_out , opt_num ) ) ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num ) )...
                                      - q_Load(k , opt_num )+ q_Solar(5 , opt_num ) == sum( x_qij( node_out , opt_num ) ) ];
                                  
     elseif( k == pv6 )
        Constraints = [ Constraints , sum( x_pij( node_in , opt_num ) )...
                                      + p_Solar( k , opt_num )-PV_cut( 6 , opt_num ) + p_Wind( k , opt_num ) - p_Load( k , opt_num ) == sum( x_pij( node_out , opt_num ) ) ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num ) )...
                                      - q_Load(k , opt_num )+ q_Solar( 6 , opt_num ) == sum( x_qij( node_out , opt_num ) ) ];                                                        
    else
        Constraints = [ Constraints , sum( x_pij( node_in , opt_num )  )...
                                      + p_Solar( k , opt_num ) + p_Wind( k , opt_num ) - p_Load( k , opt_num ) == sum( x_pij( node_out , opt_num ) ) ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num )  )...
                                      - q_Load(k , opt_num ) == sum( x_qij( node_out , opt_num ) ) ]; 
    end         
end
  
%% 支路欧姆定律
for r = 2 : 32
    Constraints = [ Constraints , x_ui_square( Branch(r,2) , opt_num ) - x_ui_square( Branch(r,3) , opt_num )  ...
                   -  ( r_ij(r) * x_pij( r , opt_num ) + x_ij(r)* x_qij( r , opt_num ) ) == 0 ];
end
             
%% 运行电压约束
Constraints = [ Constraints , x_ui_square( : , opt_num ) >= a1 ];
Constraints = [ Constraints , x_ui_square( : , opt_num ) <= a2 ];

%% 风光机组无功功率容量和最值约束
for s_n=1:q_n
%------------------------容量约束---------------------------------------
Constraints = [ Constraints ,(q_Solar(s_n , opt_num )^2 + P_pv( s_n , opt_num )^2 )<=0.9*2*( S_solar(s_n)/ ( sqrt(2) ) ) * (S_solar(s_n)/ ( sqrt(2) ) ) ];
%------------------------最值约束/功率因数约束---------------------------------------
Constraints = [ Constraints ,-Q_pv(s_n) <= q_Solar(s_n , opt_num )<= Q_pv(s_n)];
end

%------------------------sop功率因数无功约束---------------------------------------
% Constraints = [ Constraints ,-Q_sop <= Q1(: , opt_num )<= Q_sop];
% Constraints = [ Constraints ,-Q_sop <= Q2(: , opt_num )<= Q_sop];

%% 光伏机组功率削减约束
Constraints = [ Constraints , 0<= PV_cut( 1 , opt_num )<= p_Solar(pv1 , opt_num )];
Constraints = [ Constraints , 0<= PV_cut( 2 , opt_num )<= p_Solar(pv2 , opt_num )];
Constraints = [ Constraints , 0<= PV_cut( 3 , opt_num )<= p_Solar(pv3 , opt_num )];
Constraints = [ Constraints , 0<= PV_cut( 4 , opt_num )<= p_Solar(pv4 , opt_num )];
Constraints = [ Constraints , 0<= PV_cut( 5 , opt_num )<= p_Solar(pv5 , opt_num )];
Constraints = [ Constraints , 0<= PV_cut( 6 , opt_num )<= p_Solar(pv6 , opt_num )];


Constraints = [ Constraints , P_pv( 1 , opt_num )== p_Solar(pv1 , opt_num )-PV_cut( 1 , opt_num )];
Constraints = [ Constraints , P_pv( 2 , opt_num )== p_Solar(pv2 , opt_num )-PV_cut( 2 , opt_num )];
Constraints = [ Constraints , P_pv( 3 , opt_num )== p_Solar(pv3 , opt_num )-PV_cut( 3 , opt_num )];
Constraints = [ Constraints , P_pv( 4 , opt_num )== p_Solar(pv4 , opt_num )-PV_cut( 4 , opt_num )];
Constraints = [ Constraints , P_pv( 5 , opt_num )== p_Solar(pv5 , opt_num )-PV_cut( 5 , opt_num )];
Constraints = [ Constraints , P_pv( 6 , opt_num ) == p_Solar(pv6 , opt_num )-PV_cut( 6 , opt_num )];

%% 变压器1-2节点电压变比约束
% Constraints = [ Constraints , Ktij( opt_num ) == (kk - K_max) * btij( : , opt_num ) ];
% Constraints = [ Constraints , sum( btij( : , opt_num ) ) == 1 ];
Constraints = [ Constraints , x_ui_square( 2 , opt_num )  == ( k_ij0 + Ktij( opt_num ) * delta_kij ) * x_ui_square( 1 , opt_num )];

Constraints = [ Constraints , Aux( : , opt_num ) >= x_ui_square( : , opt_num )  - 1.03 ];
Constraints = [ Constraints , Aux( : , opt_num ) >= -x_ui_square( : , opt_num )  + 0.97];
Constraints = [ Constraints , Aux( : , opt_num ) >= 0 ];
toc
end

%% 储能约束
%-------------------------新加部分储能的连续性----------------
if (t==1)
 Soc_ess(1,1)=0.5;
else
 Soc_ess(1,1)=S_ESS(t-1,2);
end

%% 初始运行点ES-SOP运行约束
for opt=1:N
for rr=1:2
%ACDC损耗
Constraints = [ Constraints , P1_L( rr , opt) == C_a0*b1( rr , opt ) + C_a1*k11( rr , opt) + C_a2*k12(rr , opt) ];
Constraints = [ Constraints , P2_L( rr , opt ) == C_a0*b2( rr , opt ) + C_a1*k21( rr, opt ) + C_a2*k22( rr , opt) ];
%ACDC辅助变量定义松弛
Constraints = [ Constraints , norm( [P1( rr , opt) Q1( rr , opt )])<= k11( rr , opt) ];
Constraints = [ Constraints , norm( [2*P1( rr , opt ) 2*Q1( rr , opt) (1-k12( rr , opt))])<= 1+ k12( rr , opt)];
Constraints = [ Constraints , norm( [P2( rr , opt) Q2( rr , opt)])<= k21( rr , opt) ];
Constraints = [ Constraints , norm( [2*P2( rr , opt) 2*Q2( rr , opt ) (1-k22( rr , opt))])<= 1+ k22( rr , opt )];
%ACDC切断控制
Constraints = [ Constraints ,k11( rr , opt )<=b1( rr , opt)*S_acdc ];
Constraints = [ Constraints ,k12( rr , opt )<=b1( rr , opt )*(S_acdc)^2 ];
Constraints = [ Constraints ,k21( rr , opt)<=b2( rr , opt )*S_acdc ];
Constraints = [ Constraints ,k22( rr , opt )<=b2( rr , opt )*(S_acdc)^2 ];

%容量约束
Constraints = [ Constraints , norm( [P1( rr , opt) Q1( rr , opt )])<= 0.9*S_acdc ];
Constraints = [ Constraints , norm( [P2( rr , opt) Q2( rr , opt )])<=0.9*S_acdc ];
end

%DCDC损耗
Constraints = [ Constraints , P3_L( 1 , opt ) == C_d0*b3( 1 , opt) + C_d1*P3_abs( 1 , opt ) + C_d2*k3(1 , opt ) ];
%DCDC辅助变量定义松弛
Constraints = [ Constraints ,P3( 1 , opt)<=P3_abs( 1 , opt) ];
Constraints = [ Constraints ,-P3( 1 , opt)<=P3_abs( 1 , opt) ];
Constraints = [ Constraints , norm( [2*P3( 1 , opt) (1-k3( 1 , opt))])<= 1+ k3( 1 , opt )];
%DCDC切断控制
Constraints = [ Constraints ,P3_abs( 1 , opt )<=b3( 1 , opt )*S_dcdc ];
Constraints = [ Constraints ,k3( 1 , opt )<=b3( 1 , opt )*(S_dcdc)^2 ];

%ESS损耗
Constraints = [ Constraints , P4_L( 1 , opt ) == (C_aux + C_ES*k4( 1 , opt ))];
Constraints = [ Constraints , norm( [2*P4( 1 , opt ) (1-k4( 1 , opt ))])<= 1+ k4( 1 , opt )];

%储能运行约束
P_max=0.6;
S_min=0.2;
S_max=1;
Constraints = [ Constraints , Soc_ess( 1 , opt+1 )== Soc_ess( 1 , opt )-( P4_L( 1 , opt ) + P4( 1 , opt) )*0.25];
Constraints = [ Constraints ,-P_max <= P4( 1 , opt ) ];
Constraints = [ Constraints , P4( 1 , opt ) <= P_max ];

%容量约束
Constraints = [ Constraints ,-S_dcdc <= P3( 1 , opt ) ];
Constraints = [ Constraints , P3( 1 , opt ) <=S_dcdc ];

%有功平衡
Constraints = [ Constraints , P3( 1 , opt ) == -(P4( 1 , opt )-P4_L( 1 , opt ))];
Constraints = [ Constraints , sum(P1( : , opt ) + P2( : , opt)+ P1_L( : , opt )+ P2_L( : , opt ))+ P3( 1 , opt ) + P3_L( 1 , opt )==0];


%多馈线约束
Constraints = [ Constraints , sum(b1( : , opt )) <= 1 ];
Constraints = [ Constraints , sum(b1( 1 , opt )+ b2( 2 , opt )) <= 1 ];
Constraints = [ Constraints , sum(b2( : , opt ))  <= 1 ];
Constraints = [ Constraints , sum(b1( 2 , opt )+ b2( 1 , opt ))<= 1 ];
end

Constraints = [ Constraints ,S_min <= Soc_ess( 1 , : )  ];
Constraints = [ Constraints ,Soc_ess( 1 , : ) <= S_max ];

%----------------目标函数求和------------------------------------------
f1 = sum(f);

%----------------求解器及求解设置------------------------------------------
options = sdpsettings('verbose',1,'solver','gurobi');%gurobi
sol = solvesdp(Constraints,f1,options);
% --------------------------Analyze error flags----------------------------
if (sol.problem == 0)
 % Extract and display value
    result = double(f1)
else
    
 display('Hmm, something went wrong!');
 sol.info
 yalmiperror(sol.problem)
end
S_ESS(t,:)= Soc_ess(1,:);%每次求解的荷电状态
%% 数据保存
%------------------------电压、电流、支路有功无功--------------------------
x_ui_squaret(:,t)=x_ui_square(:,1);
x_pijt(:,t)=x_pij(:,1); 
x_qijt(:,t)=x_qij(:,1); 
%-----------------------损耗-----------------------------------
f_startt(:,t)=f_start(:,1); 
f_essopt(:,t)=f_essop(:,1); 
f_PV_cutt(:,t)=f_PV_cut(:,1); 
%-----------------------ES-SOP-----------------------------------
b1t(:,t)=b1(:,1);  
b2t(:,t)=b2(:,1);   
b3t(:,t)=b3(:,1);    
%AC/DC视在功率辅助变量/非负
k11t(:,t)=k11(:,1); 
k12t(:,t)=k12(:,1); 
k21t(:,t)=k21(:,1); 
k22t(:,t)=k22(:,1); 
%四个子系统损耗
P1_Lt(:,t)=P1_L(:,1); 
P2_Lt(:,t)=P2_L(:,1); 
P3_Lt(:,t)=P3_L(:,1); 
P4_Lt(:,t)=P4_L(:,1); 
%四个子系统有功功率
P1t(:,t)=P1(:,1); 
P2t(:,t)=P2(:,1); 
P3t(:,t)=P3(:,1); 
P4t(:,t)=P4(:,1); 
 
%AC/DC无功功率
Q1t(:,t)=Q1(:,1); 
Q2t(:,t)=Q2(:,1); 
%DC/DC损耗中的绝对值项和平方项
P3_abst(:,t)=P3_abs(:,1); 
k3t(:,t)=k3(:,1);
%ESS损耗中的平方项和荷电状态
k4t(:,t)=k4(:,1); 
Soc_esst(:,t)=Soc_ess(:,2);
%-----------------------光伏无功输出-----------------------------------
q_Solart(:,t)=q_Solar(:,1);
PV_cutt(:,t)=PV_cut(:,1);
%-----------------------保存96组sop有功无功、pv无功-----------------------------------
b1_1(t,:)=b1(1,:); 
b1_2(t,:)=b1(2,:); 
b2_1(t,:)=b2(1,:); 
b2_2(t,:)=b2(2,:); 
b3_1(t,:)=b3(1,:); 

P1_12(t,:)=P1(1,:); 
P1_25(t,:)=P1(2,:); 
P2_22(t,:)=P2(1,:); 
P2_29(t,:)=P2(2,:); 

Q1_12(t,:)=Q1(1,:); 
Q1_25(t,:)=Q1(2,:); 
Q2_22(t,:)=Q2(1,:); 
Q2_29(t,:)=Q2(2,:); 

q_Solar_10(t,:)=q_Solar(1,:); 
q_Solar_24(t,:)=q_Solar(2,:); 
q_Solar_28(t,:)=q_Solar(3,:); 
% q_Solar_24(t,:)=q_Solar(4,:); 
% q_Solar_28(t,:)=q_Solar(5,:); 
end
disp(['求解完成',num2str(etime(clock,t1))]);

%% 数据读取
%------------------------电压、电流、支路有功无功--------------------------
x_ui_square_data20 = double((x_ui_squaret));
x_pij_data20 = double(x_pijt);
x_qij_data20 = double(x_qijt);
%-----------------------损耗-----------------------------------
f_start_data20=double(f_startt);
f_essop_data20=double(f_essopt);
f_PV_cut_data20=double(f_PV_cutt);
%-----------------------ES-SOP-----------------------------------
%VSC接通/断开变量
b1_data20= double(b1t);
b2_data20= double(b2t);
b3_data20= double(b3t);
%AC/DC视在功率辅助变量/非负
k11_data20= double(k11t);
k12_data20= double(k12t);
k21_data20= double(k21t);
k22_data20= double(k22t);
%四个子系统损耗
P1_L_data20=double(P1_Lt);
P2_L_data20=double(P2_Lt);
P3_L_data20=double(P3_Lt);
P4_L_data20=double(P4_Lt);
%四个子系统有功功率
P1_data20=double(P1t);
P2_data20=double(P2t);
P3_data20=double(P3t);
P4_data20=double(P4t);
%AC/DC无功功率
Q1_data20=double(Q1t);
Q2_data20=double(Q2t);
%DC/DC损耗中的绝对值项和平方项
P3_abs_data20=double(P3_abst);
k3_data20=double(k3t);
%ESS损耗中的平方项和荷电状态
k4_data20=double(k4t);  
Soc_ess_data20=double(Soc_esst);
%-----------------------电压偏差-----------------------------------
aux20=sum(sum(abs(x_ui_square_data20-1)))*0.25;
%-----------------------光伏无功输出-----------------------------------
q_Solar_data20=double(q_Solart);
PV_cut_data20=double(PV_cutt);
%-----------------------保存96组sop有功无功、pv无功-----------------------------------
b1_1_9620=double(b1_1);
b1_2_9620=double(b1_2);
b2_1_9620=double(b2_1);
b2_2_9620=double(b2_2);
b3_1_9620=double(b3_1);

P1_12_9620=double(P1_12);
P1_25_9620=double(P1_25);
P2_22_9620=double(P2_22);
P2_29_9620=double(P2_29);

Q1_12_9620=double(Q1_12);
Q1_25_9620=double(Q1_25);
Q2_22_9620=double(Q2_22);
Q2_29_9620=double(Q2_29);

q_Solar_10_9620=double(q_Solar_10);
% q_Solar_11_9620=double(q_Solar_11);
% q_Solar_21_9620=double(q_Solar_21);
q_Solar_24_9620=double(q_Solar_24);
q_Solar_28_9620=double(q_Solar_28);

%% 存储数据
save stage2_1_WithoutVSM_随机性 x_ui_square_data20  x_pij_data20 x_qij_data20...
f_start_data20  f_essop_data20  f_PV_cut_data20...
k11_data20 k21_data20...
P1_data20 P2_data20 P3_data20 P4_data20 ...
P1_L_data20 P2_L_data20 P3_L_data20 P4_L_data20...
Q1_data20 Q2_data20 ...
Soc_ess_data20 ...
aux20 ...
q_Solar_data20 PV_cut_data20...
P1_12_9620 P1_25_9620 P2_22_9620 P2_29_9620...
Q1_12_9620 Q1_25_9620 Q2_22_9620 Q2_29_9620...
q_Solar_10_9620  q_Solar_24_9620  q_Solar_28_9620...
b1_1_9620 b1_2_9620 b2_1_9620 b2_2_9620 b3_1_9620...
b1_data20 b2_data20 b3_data20;


% save 日内步骤1传递给2的数据1_96组_WithoutVSM_随机性 P1_12_9620 P1_25_9620 P2_22_9620 P2_29_9620...
% Q1_12_9620 Q1_25_9620 Q2_22_9620 Q2_29_9620...
% q_Solar_10_9620  q_Solar_24_9620 q_Solar_28_9620 ...
% b1_1_9620 b1_2_9620 b2_1_9620 b2_2_9620 b3_1_9620;

a1=sum(P3_L_data20+P4_L_data20)*0.25
a2=sum(f_essop_data20)





