%% SDM283 - MiniProject1 - IRB460
%% Task1 - by Group8 MIAO Ziliang 11911901
clc; clear all; close all;
theta1 = 25;
theta2 = 40;
theta3 = 20;
global q00 q10 q20 q30 gab

%输入各个joint
q00 = 0.001*[0; 0; 0];
q10 = 0.001*[0; 260; 742.5];
q20 = 0.001*[0; 260; 742.5+945];
q30 = 0.001*[0; 260+1025; 742.5+945];

%额外可能用到的joint
q40 = 0.001*[0; 260+1025+220; 742.5+945-217.5];
q50 = 0.001*[0; -140; 742.5];
q60 = 0.001*[0; 0; 742.5+140.84];
q70 = 0.001*[0; 260+229.81; 742.5+945+192.84];
q80 = 0.001*[0; 0; 742.5+140.84+945];
q90 = 0.001*[0; 260+229.81+1025; 742.5+945+192.84];
q100 = 0.001*[0; -140; 742.5+945];

%各个杆的重心和三角形重心
c10 = [0;0.1479760992;0.5686888371];
c20 = [0;0.2834415217;1.1379921239];
c30 = [0;0.136553367;0.7422922859];
c40 = [0;-0.1400000011;1.2150000115];
c50 = [0;0.4379407461;1.7370254708];
c60 = [0;-0.034408645058;1.3558415391];
c70 = [0;0.2535938282;1.770122769];
c80 = [0;1.0023172834;1.8803362521];
c90 = [0;1.4437790926;1.6191476398];

gab = CoordinateTransformation(theta1,theta2,theta3);

%以下操作为了方便ADAMS仿真验证（提取变换矩阵的最后一列中的前三行：坐标）
%通过变换矩阵对end effector的初始坐标进行变换
Q4 = gab * [0;0;0;1];
q4 = Q4(1:3,:);

%输出变换矩阵
disp(gab);
fprintf("Coordinate of the end of IRB460 is [%.4f,%.4f,%.4f] at the condition of theta1 = %.4f, theta2 = %.4f, theta3 = %.4f\n",q40(1),q40(2),q40(3),theta1,theta2,theta3);

%% Task2 - by Group8 MIAO Ziliang 11911901
%输入各关节角度变化速率
theta1_dot = 0;
theta2_dot = 8;
theta3_dot = 4;

%在此情况下计算
Vab_hat1 = SpatialVelocity(theta1, theta2, theta3, theta1_dot, theta2_dot, theta3_dot);

%以下操作为了方便ADAMS仿真验证（提取变换矩阵的最后一列中的前三行：坐标）
V4 = Vab_hat1 * [q40;1];
v4 = V4(1:3,:);

%% Task3 - by Group8 ZHOU JinHui 11911908
%输入payload
F = [0;0;-800;0;0;0];

 %输入节点旋转轴
 x = [1; 0; 0]; 
 z = [0; 0; 1];
    
%输入单位角速度
x_hat = [0 0 0; 0 0 -1; 0 1 0];
z_hat = [0 -1 0; 1 0 0; 0 0 0];
 
%计算旋量
twist1 = [z_hat, cross(q00,z);0 0 0 0];
twist2 = [x_hat, cross(q10,x);0 0 0 0];
twist3 = [x_hat, cross(q20,x);0 0 0 0];
twist4 = [x_hat, cross(q30,x);0 0 0 0];
    
%计算伴随矩阵
Ad1 = Adjoint(RotationZ(q00, theta1));
Ad2 = Adjoint(RotationX(q10,theta2));
Ad3 = Adjoint(RotationX(q20,(theta3-theta2)));

twist1_ = [cross(q00,z);z];
twist2_ = [cross(q10,x);x];
twist3_ = [cross(q20,x);x];
twist4_ = [cross(q30,x);x];

%计算“Twist in new location”
twist1_new = twist1_;
twist2_new = (Ad1) * twist2_;
twist3_new = (Ad1*Ad2) * twist3_;
twist4_new = (Ad1*Ad2*Ad3) * twist4_;
    
%输出雅可比矩阵(需要简化为三个自由度——末端为旋转)
%J = [twist1_new,twist2_new,twist3_new,twist4_new];
J = [twist1_new,twist2_new-twist3_new,twist3_new-twist4_new];
%此处可能需要对矩阵进行化简（至三个自由度）
disp(J);
Vab1 = J*[theta1_dot;theta2_dot;theta3_dot];
Vab_hat2=VtoVha(Vab1);
% Vab_hat2 =[[0,-Vab(6),Vab(5);Vab(6),0,-Vab(4);-Vab(5),Vab(4),0],[Vab(1);Vab(2);Vab(3)];0,0,0,0];
V42 = Vab_hat2 * [q40;1];


%计算各个质心的坐标表达式
syms theta1 theta2 theta3
c1 = RotationZ(c10,theta1)*[c10;1];
c2 = RotationZ(c10,theta1)*RotationX(q10,theta2)*[c20;1];
c3 = RotationZ(c10,theta1)*RotationX(q10,theta3)*[c30;1];
c4 = RotationZ(c10,theta1)*RotationX(q10,theta3)*RotationX(q50,theta2-theta3)*[c40;1];
c5 = RotationZ(c10,theta1)*RotationX(q10,theta2)*RotationX(q20,-theta2+theta3)*[c50;1];
c6 = RotationZ(c10,theta1)*RotationX(q60,theta2)*[c60;1];
c7 = RotationZ(c10,theta1)*RotationX(q10,theta2)*RotationX(q20,-theta2)*[c70;1];
c8 = RotationZ(c10,theta1)*RotationX(q10,theta2)*RotationX(q20,-theta2)*RotationX(q70,theta3)*[c80;1];
c9 = RotationZ(c10,theta1)*RotationX(q10,theta2)*RotationX(q20,-theta2+theta3)*RotationX(q30,-theta3)*[c90;1];

%计算静态下重力势能
%先声明各个质量
m1 = 102.7282598605;
m2 = 12.4193410149;
m3 = 5.2824118449;
m4 = 2.9463248039;
m5 = 9.6977245306;
m6 = 1.4433507166;
m7 = 2.8514044732;
m8 = 1.4155689868;
m9 = 6.0134646985;
G = 9.80665;
V = G*(m1*[0;0;1;0]'*c1 + m2*[0;0;1;0]'*c2 + m3*[0;0;1;0]'*c3 + m4*[0;0;1;0]'*c4 + m5*[0;0;1;0]'*c5 + ...
    + m6*[0;0;1;0]'*c6 + m7*[0;0;1;0]'*c7 + m8*[0;0;1;0]'*c8 + m9*[0;0;1;0]'*c9);

%计算力矩
V_diff_theta1 = diff(V,theta1);
V_diff_theta2 = diff(V,theta2);
V_diff_theta3 = diff(V,theta3);
V_diff_theta = [V_diff_theta1;V_diff_theta2;V_diff_theta3];
torque_nogravity = J'*F;
torque_gravity = -V_diff_theta + J'*F;
%赋值
vpa(subs(torque_gravity,{theta1,theta2,theta3},{0,0,0}))

%Task4 - by Group8 Liu HongLei 11912905

%计算动能
%先声明各个joint角速度
syms dot_theta1 dot_theta2 dot_theta3
dot_c1 = diff(c1,theta1)*dot_theta1 + diff(c1,theta2)*dot_theta2 + diff(c1,theta3)*dot_theta3;
dot_c2 = diff(c2,theta1)*dot_theta1 + diff(c2,theta2)*dot_theta2 + diff(c2,theta3)*dot_theta3;
dot_c3 = diff(c3,theta1)*dot_theta1 + diff(c3,theta2)*dot_theta2 + diff(c3,theta3)*dot_theta3;
dot_c4 = diff(c4,theta1)*dot_theta1 + diff(c4,theta2)*dot_theta2 + diff(c4,theta3)*dot_theta3;
dot_c5 = diff(c5,theta1)*dot_theta1 + diff(c5,theta2)*dot_theta2 + diff(c5,theta3)*dot_theta3;
dot_c6 = diff(c6,theta1)*dot_theta1 + diff(c6,theta2)*dot_theta2 + diff(c6,theta3)*dot_theta3;
dot_c7 = diff(c7,theta1)*dot_theta1 + diff(c7,theta2)*dot_theta2 + diff(c7,theta3)*dot_theta3;
dot_c8 = diff(c8,theta1)*dot_theta1 + diff(c8,theta2)*dot_theta2 + diff(c8,theta3)*dot_theta3;
dot_c9 = diff(c9,theta1)*dot_theta1 + diff(c9,theta2)*dot_theta2 + diff(c9,theta3)*dot_theta3;

T = 0.5*m1*(dot_c1)'*(dot_c1) + 0.5*m2*(dot_c2)'*(dot_c2) + 0.5*m3*(dot_c3)'*(dot_c3) + 0.5*m4*(dot_c4)'*(dot_c4) + 0.5*m5*(dot_c5)'*(dot_c5) + 0.5*m6*(dot_c6)'*(dot_c6) + 0.5*m7*(dot_c7)'*(dot_c7) + 0.5*m8*(dot_c8)'*(dot_c8) + 0.5*m9*(dot_c9)'*(dot_c9);
%动能对角速度求导
T_diff_dot_theta1 = diff(T,dot_theta1);
T_diff_dot_theta2 = diff(T,dot_theta2);
T_diff_dot_theta3 = diff(T,dot_theta3);
T_diff_dot_theta = [T_diff_dot_theta1;T_diff_dot_theta2;T_diff_dot_theta3;];
%动能对角速度求导再对t求导
%先声明角加速度
ddot_theta1 = 0*pi/180;
ddot_theta2 = 1*pi/180;
ddot_theta3 = 0*pi/180;
T_diff_t1 = diff(T_diff_dot_theta1,theta1)*dot_theta1 + diff(T_diff_dot_theta1,theta2)*dot_theta2 + ...
    + diff(T_diff_dot_theta1,theta3)*dot_theta3 + diff(T_diff_dot_theta1,dot_theta1)*ddot_theta1 + ... 
    + diff(T_diff_dot_theta1,dot_theta2)*ddot_theta2 + diff(T_diff_dot_theta1,dot_theta3)*ddot_theta3;
T_diff_t2 = diff(T_diff_dot_theta2,theta1)*dot_theta1 + diff(T_diff_dot_theta2,theta2)*dot_theta2 + ...
    + diff(T_diff_dot_theta2,theta3)*dot_theta3 + diff(T_diff_dot_theta2,dot_theta1)*ddot_theta1 + ... 
    + diff(T_diff_dot_theta2,dot_theta2)*ddot_theta2 + diff(T_diff_dot_theta2,dot_theta3)*ddot_theta3;
T_diff_t3 = diff(T_diff_dot_theta3,theta1)*dot_theta1 + diff(T_diff_dot_theta3,theta2)*dot_theta2 + ...
    + diff(T_diff_dot_theta3,theta3)*dot_theta3 + diff(T_diff_dot_theta3,dot_theta1)*ddot_theta1 + ... 
    + diff(T_diff_dot_theta3,dot_theta2)*ddot_theta2 + diff(T_diff_dot_theta3,dot_theta3)*ddot_theta3;
T_diff_t = [T_diff_t1;T_diff_t2;T_diff_t3];

%动能对角度求导
T_diff_theta1 = diff(T,theta1);
T_diff_theta2 = diff(T,theta2);
T_diff_theta3 = diff(T,theta3);
T_diff_theta = [T_diff_theta1;T_diff_theta2;T_diff_theta3;];
%力矩
torque_dynamics = T_diff_t - T_diff_theta + V_diff_theta;
%赋值
vpa(subs(torque_dynamics,{theta1,theta2,theta3,dot_theta1,dot_theta2,dot_theta3},{0*pi/180,17.5*pi/180,0*pi/180,0*pi/180,1*pi/180,0*pi/180}))


%空间速度计算公式
function Vab = SpatialVelocity(Theta1, Theta2, Theta3, Theta1_dot, Theta2_dot, Theta3_dot)
    %声明变量
    global q00 q10 q20 q30 
    %进行角度弧度制转换
    Theta1 = Theta1/180*pi;
    Theta2 = Theta2/180*pi;
    Theta3 = Theta3/180*pi; 
    Theta1_dot = Theta1_dot/180*pi;
    Theta2_dot = Theta2_dot/180*pi;
    Theta3_dot = Theta3_dot/180*pi;
    
    %输入节点旋转轴
    x = [1; 0; 0]; 
    z = [0; 0; 1];
    
    %输入单位角速度
    x_hat = [0 0 0; 0 0 -1; 0 1 0];
    z_hat = [0 -1 0; 1 0 0; 0 0 0];
    
    %计算旋转矩阵
    R1 = RotationZ(q00, Theta1);
    R2 = RotationX(q10,Theta2);
    R3 = RotationX(q20,(Theta3-Theta2));
    R4 = RotationX(q30,(-Theta3));
    
    %计算旋量
    twist1 = [z_hat, cross(q00,z);0 0 0 0];
    twist2 = [x_hat, cross(q10,x);0 0 0 0];
    twist3 = [x_hat, cross(q20,x);0 0 0 0];
    twist4 = [x_hat, cross(q30,x);0 0 0 0];
    
    %计算gab0、gab_dot、gab_inverse
    gab0 = [eye(3), [0;1025+260+220;945+742.5-217.5];0 0 0 1];
    gab_dot = (twist1 * Theta1_dot * R1 * R2 * R3 * R4... 
    + R1 * twist2 * Theta2_dot * R2 * R3 * R4... 
    +R1 * R2 * twist3 * (Theta3_dot - Theta2_dot) * R3 * R4... 
    + R1 * R2 * R3 * twist4 * (-Theta3_dot) * R4) * gab0;
    gab_inverse = inv(gab0) * inv(R4) * inv(R3) * inv(R2) * inv(R1);
    
    %最后计算空间速度
%     Vab = gab_dot * gab_inverse;
    Vab = twist1 * Theta1_dot + R1*twist2*Theta2_dot*inv(R1) + ...
        R1*R2*twist3*(Theta3_dot-Theta2_dot)*inv(R2)*inv(R1) + ...
        R1*R2*R3*twist4*(-Theta3_dot)*inv(R3)*inv(R2)*inv(R1);
end

%坐标系变换函数
function T = CoordinateTransformation(theta1, theta2, theta3)
    global q00 q10 q20 q30
    
    %进行角度弧度制转换
    theta1 = theta1/180*pi;
    theta2 = theta2/180*pi;
    theta3 = theta3/180*pi;
    
    %输入end effector的初始坐标
    gab0 = [eye(3), 0.001*[0;1025+260+220;945+742.5-217.5];0 0 0 1];
    
    %引用Rotation函数进行变换叠加
    gab = RotationZ(q00, theta1) * RotationX(q10,theta2) * RotationX(q20,(theta3-theta2)) * RotationX(q30,(-theta3)) * gab0;
    T = gab;
end

%计算twist:X
% function t = twist(q_axis)
%  %x代表旋转轴
%     x = [1; 0; 0];
%     %x_hat代表单位角速度
%     x_hat = [0 0 0; 0 0 -1; 0 1 0];
%     
%     %%Caculate Method 1 (Use method "expm")：利用指数函数进行计算
%     t = [x_hat, cross(q_axis,x); 0 0 0 0];
% end

%旋转矩阵计算函数：X
function R = RotationX(q_axis,theta)
    %x代表旋转轴
    x = [1; 0; 0];
    %x_hat代表单位角速度
    x_hat = [0 0 0; 0 0 -1; 0 1 0];
    
    %%Caculate Method 1 (Use method "expm")：利用指数函数进行计算
    %twist = [x_hat, cross(q_axis,x); 0 0 0 0];
    %R = expm(twist * theta);
    
    %%Caculate Method 2 (Simplify by Rodrigues Formula)：直接利用Rodriguez公式计算
    r = eye(3) + x_hat * sin(theta) + (x_hat)^2 * (1 - cos(theta));
    R = [r, (eye(3) - r) * q_axis; 0 0 0 1];
end

%旋转矩阵计算函数：Z
function R = RotationZ(q_axis,theta)
    z = [0; 0; 1];
    z_hat = [0 -1 0; 1 0 0; 0 0 0];
    %Using Method 2 here
    r = eye(3) + z_hat * sin(theta) + (z_hat)^2 * (1 - cos(theta));
    R = [r, (eye(3) - r) * q_axis; 0 0 0 1];
end

%计算对应的伴随矩阵
function A = Adjoint(Rotation)
    R = Rotation([1,2,3],[1,2,3]);
    pI = Rotation([1,2,3],(4));
    p = [0,-pI(3),pI(2);pI(3),0,-pI(1);-pI(2),pI(1),0];
    A = [R,p*R;zeros(3),R];
end

function Vhat = VtoVhat(V)
    v = [V(1,1);V(2,1);V(3,1)];
    w = [V(4,1);V(5,1);V(6,1)];
    Vhat = [vectorTohat(w) v; 0 0 0 0];
end
);


