function [B] = get_inertia_matrix(q, params)
% Funzione per ricavare la matrice di inerzia B.

% q: vettore delle coordinate lagrangiane (theta1 theta2)
% params: vettore dei parametri caratteristici del robot (m1 a1 l1 I1 m2 a2 l2 I2)

% Definisco variabili di interesse:
m1 = params(1);
a1 = params(2);
l1 = params(3);
I1 = params(4);
m2 = params(5);
a2 = params(6);
l2 = params(7);
I2 = params(8);


%Jacobiani necessari per il calcolo della matrice d'inerzia:
Jp1 = [-l1.*sin(q(1)), 0;l1.*cos(q(1)), 0; 0, 0];
Jp2 = [-a1.*sin(q(1))-l2.*sin(q(1)+q(2)), -l2.*sin(q(1)+q(2)); ...
    a1.*cos(q(1))+l2.*cos(q(1)+q(2)), l2.*cos(q(1)+q(2)); 0, 0];
J01 = [0, 0; 0, 0; 1, 0];
J02 = [0, 0; 0, 0; 1, 1];


%Calcolo della matrice d'inerzia:
B = m1*transpose(Jp1)*Jp1 + m2*transpose(Jp2)*Jp2 ...
    + I1*transpose(J01)*J01 + I2*transpose(J02)*J02;
end