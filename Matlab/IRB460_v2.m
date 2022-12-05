clc; clear all; close all;
syms theta1 theta2 theta3
p1 = [0; 0; 234.5];
p2 = [0; 260; 508];
p3 = [0; 0; 0];
p4 = [0; 945;0];
p5 = [1025; 0; 0];

T1 = [RotationZ(theta1),p1;0 0 0 1];
T2 = [RotationY(pi/2),p2;0 0 0 1];
T3 = [RotationZ(pi/2-theta2+theta3),p3;0 0 0 1];
T4 = [RotationZ(pi),p4;0 0 0 1];
T5 = [RotationZ(135.33/180*pi-theta3),p5;0 0 0 1];
T = T1 * T2 * T3 * T4 * T5;
subs(T,[theta1,theta2,theta3],[0,0,0])

exp([1;0;0]*theta2)

function R = RotationZ(theta)
    R = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
end
function R = RotationY(theta)
    R = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
end