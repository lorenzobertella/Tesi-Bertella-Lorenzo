function [B] = get_inertia_matrix(q, params)
%Funzione per ricavare la matrice di inerzia B.

%q: vettore delle coordinate lagrangiane (theta1 theta2)
%params: vettore dei parametri caratteristici del robot (m1 a1 l1 I1 m2 a2 l2 I2)

%Jacobiani necessari per il calcolo della matrice d'inerzia:
Jp1 = [-params(3).*sin(q(1)) 0;params(3).*cos(q(1)) 0; 0 0];
Jp2 = [-params(2).*sin(q(1))-params(7).*sin(q(1)+q(2)) -params(7).*sin(q(1)+q(2)); params(2).*cos(q(1))+params(7).*cos(q(1)+q(2)) params(7).*cos(q(1)+q(2)); 0 0];
J01 = [0 0; 0 0; 1 0];
J02 = [0 0; 0 0; 1 1];


%Calcolo della matrice d'inerzia:
B = params(1)*transpose(Jp1)*Jp1 + params(5)*transpose(Jp2)*Jp2 ...
    + params(4)*transpose(J01)*J01 + params(8)*transpose(J02)*J02;
end