clc; clear all; close all;
q10 = [sqrt(3);1;0];
v1 = [[0 -1 0; 1 0 0; 0 0 0 ] cross(q10,[0;0;1]); 0 0 0 0];
V1 = expm(v1*pi/6);
gab0 = [eye(3), [2+sqrt(3);1;0]; 0 0 0 1];
gab = V1 * gab0;
gab * [0;0;0;1]
2*sqrt(3)