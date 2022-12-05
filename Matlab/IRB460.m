clc; clear all;
syms theta1 theta2 theta3
p1 = [0; 0; 234.5];
p2 = [0; 260; 508];
p3 = [0; 0; 0];
p4 = [0; 945;0];
p5 = [1025; 0; 0];


T_01=transformation(theta1,'z',p1);
T_12=transformation(pi/2,'y',p2);
T_23=transformation(pi/2-theta2+theta3,'z',p3);
T_34=transformation(pi,'z',p4);
T_45=transformation(135.33/180*pi-theta3,'z',p5);
T_05 = T_01 * T_12 * T_23 * T_34 * T_45;

theta1 = 0;
theta2 = 0;
theta3 = 0;
T = subs(T_05)
format short

function R = rotation(theta,direction)
    if(direction == 'x')
        R = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
    end
    if(direction == 'y')
        R = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
    end
    if(direction == 'z')
        R = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
    end
    if(direction == 's')
        R = [1 0 0; 0 1 0; 0 0 1];
    end
end
function T = transformation(theta,direction,p)
    R = rotation(theta,direction);
    T = [R p;0 0 0 1];
end
