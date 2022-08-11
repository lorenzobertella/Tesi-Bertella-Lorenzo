function [dx_dt] = ode_func(t,x,torque_control, B, C, g)
%equazioni differenziali scritte in forma di stato

%definisco le variabili di stato:
x1 = x(1);              %theta1
x2 = x(2);              %theta2
x3 = x(3);              %dtheta1_dt
x4 = x(4);              %dtheta2_dt

%calcolo intermedio
z = B\(torque_control-C*[x3;x4]-g);

%scrittura in forma di stato:
dx1_dt = x3;
dx2_dt = x4;
dx3_dt = z(1,:);
dx4_dt = z(2,:);

%uscita della ode_func:
dx_dt = [dx1_dt; dx2_dt; dx3_dt; dx4_dt];

%questa non serve più perchè ho cambiato 
%metodo di integrazione
end