function [Pe] = direct_kinematics(q,params)
%Funzione che restituisce la posizione dell'end-effector dato il vettore
%di coordinate lagrangiane

%Definizione di variabili comode:
theta1 = q(1);
theta2 = q(2);
a1 = params(2);
a2 = params(6);

%Calcolo coordinate end-effector:
Pe = [a1 .* cos(theta1) + a2 .* cos(theta1+theta2); 
    a1 .* sin(theta1) + a2 .* sin(theta1+theta2)];
end