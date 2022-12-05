%% SDM283 - MiniProject1 - IRB460
%% Taks1 - by Group8 MIAO Ziliang 11911901
clc; clear all; close all;
theta1 = 0;
theta2 = 40;
theta3 = 0;
global q00 q10 q20 q30 V1 V2
q00 = [0; 0; 0];
q10 = [0; 260; 742.5];
q20 = [0; 260; 742.5+945];
q30 = [0; 260+1025; 742.5+945]; 

Q4 = CoordinateTransformation(theta1,theta2,theta3) * [0;0;0;1];
q4 = Q4(1:3,:);
fprintf("Coordinate of the end of IRB460 is [%.4f,%.4f,%.4f] at the condition of theta1 = %.4f, theta2 = %.4f, theta3 = %.4f\n",q4(1),q4(2),q4(3),theta1,theta2,theta3);

%% Taks2 - by Group8 MIAO Ziliang 11911901
theta1_dot = 0;
theta2_dot = 8;
theta3_dot = 0;
V = SpatialVelocity(theta1, theta2, theta3, theta1_dot, theta2_dot, theta3_dot);
V40 = V * [q4;0];
v40 = V40(1:3,:);
function V = SpatialVelocity(theta1, theta2, theta3, theta1_dot, theta2_dot, theta3_dot)
    global q00 q10 q20 q30 V1 V2
    theta1 = theta1/180*pi;
    theta2 = theta2/180*pi;
    theta3 = theta3/180*pi; 
    theta1_dot = theta1_dot/180*pi;
    theta2_dot = theta2_dot/180*pi;
    theta3_dot = theta3_dot/180*pi;
    x = [1; 0; 0];
    z = [0; 0; 1];
    R1 = RotationZ(q00, theta1);
    R2 = RotationX(q10,theta2);
    R3 = RotationX(q20,(theta3-theta2));
    R4 = RotationX(q30,(-theta3));
    x_hat = [0 0 0; 0 0 -1; 0 1 0];
    z_hat = [0 -1 0; 1 0 0; 0 0 0];
    twist1 = [z_hat, cross(q00,z);0 0 0 0];
    twist2 = [x_hat, cross(q10,x);0 0 0 0];
    twist3 = [x_hat, cross(q20,x);0 0 0 0];
    twist4 = [x_hat, cross(q30,x);0 0 0 0];
    
    gab0 = [eye(3), [0;1025+260+220;945+742.5-217.5];0 0 0 1];
    gab_dot = (twist1 * theta1_dot * R1 * R2 * R3 * R4... 
    + R1 * twist2 * theta2_dot * R2 * R3 * R4... 
    +R1 * R2 * twist3 * (theta3_dot - theta2_dot) * R3 * R4... 
    + R1 * R2 * R3 * twist4 * (-theta3_dot) * R4) * gab0;
    gab_inverse = inv(gab0) * inv(R4) * inv(R3) * inv(R2) * inv(R1);
    V1 = gab_dot * gab_inverse;
    V2 = twist1 * theta1_dot + R1*twist2*theta2_dot*inv(R1) + ...
        R1*R2*twist3*(theta3_dot-theta2_dot)*inv(R2)*inv(R1) + ...
        R1*R2*R3*twist4*(-theta3_dot)*inv(R3)*inv(R2)*inv(R1);
    V = twist2 * theta2_dot + R2 * twist3 * inv(R2) * (theta3_dot - theta2_dot) +...
        (R2 * R3)*twist4*(inv(R3)*inv(R2))*(-theta3_dot);
end


function T = CoordinateTransformation(theta1, theta2, theta3)
    global q00 q10 q20 q30
    theta1 = theta1/180*pi;
    theta2 = theta2/180*pi;
    theta3 = theta3/180*pi;  
    gab0 = [eye(3), [0;1025+260+220;945+742.5-217.5];0 0 0 1];
    gab = RotationZ(q00, theta1) * RotationX(q10,theta2) * RotationX(q20,(theta3-theta2)) * RotationX(q30,(-theta3)) * gab0;
    T = gab;
end
function R = RotationX(q_axis,theta)
    x_hat = [0 0 0; 0 0 -1; 0 1 0];
    x = [1; 0; 0];
    % Caculate Method 1 (Use method "expm")
%     twist = [x_hat, cross(q_axis,x); 0 0 0 0];
%     R = expm(twist * theta);
    % Caculate Method 2 (Simplify by Rodrigues Formula)
    r = eye(3) + x_hat * sin(theta) + (x_hat)^2 * (1 - cos(theta));
    R = [r, (eye(3) - r) * q_axis; 0 0 0 1];
end
function R = RotationZ(q_axis,theta)
    z_hat = [0 -1 0; 1 0 0; 0 0 0];
    z = [0; 0; 1];
    r = eye(3) + z_hat * sin(theta) + (z_hat)^2 * (1 - cos(theta));
    R = [r, (eye(3) - r) * q_axis; 0 0 0 1];
end


