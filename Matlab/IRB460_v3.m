
% syms theta1 theta2 theta3
theta1 = 10/180 * pi;
theta2 = 10/180*pi;
theta3 = 10/180*pi;

q10 = [0; 260; 742.5];
q20 = [0; 260; 742.5+945];
q30 = [0; 260+1025; 742.5+945];

x_hat = [0 0 0; 0 0 -1; 0 1 0];
x = [1; 0; 0];
twist1 = [x_hat, cross(q10,x); 0 0 0 0];
twist2 = [x_hat, cross(q20,x); 0 0 0 0];
twist3 = [x_hat, cross(q30,x); 0 0 0 0];

R1 = expm(twist1*theta2);
R2 = expm(twist2*(theta3 - theta2));
R3 = expm(twist3*(-theta3));


%% 验证q2坐标
g20 = [eye(3) [0; 260; 1687.5]; 0 0 0 1];
g2 = R1 * g20 ;
Q2 = g2 * [0; 0; 0; 1];
%% 验证q3坐标
g30 = [eye(3) [0;1025+260;1687.5]; 0 0 0 1];
g3 = R1 * R2 * g30;
Q3 = g3 * [0;0;0;1];
%% 计算末端坐标
gab0 = [eye(3), [0;1025+260+220;945+742.5-217.5];0 0 0 1];
gab = R1 * R2 * R3 * gab0;
Gab = [RotateZ(theta1),[0;0;0];0 0 0 1] * gab;
Q41 = Gab * [0; 0; 0; 1];


%% 直接计算坐标
 p1 = [0; 260; 742.5];
 p2 = p1 + RotateX(theta2) * [0; 0; 945];
 p3 = p2 + RotateX(theta3) * [0; 1025; 0];
 p4 = p3 + [0; 220; -217.5];



function R = RotateX(theta)
    R = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
end
function R = RotateZ(theta)
     R = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
end


