function [g] = get_gravity_vector(q, params)
%Funzione per ricavare il vettore dei termini gravitazionali g.

%q: vettore delle coordinate lagrangiane (theta1 theta2)
%params: vettore dei parametri caratteristici del robot (m1 a1 l1 I1 m2 a2 l2 I2)

%Vettore accelerazione di gravità:
g0 = [0; -9.81; 0];           %quante cifre significative???

%Jacobiani necessari per il calcolo della matrice di gravità:
Jp1 = [-params(3).*sin(q(1)) 0; params(3).*cos(q(1)) 0; 0 0];
Jp2 = [-params(2).*sin(q(1))-params(7).*sin(q(1)+q(2)) -params(7).*sin(q(1)+q(2)); params(2).*cos(q(1))+params(7).*cos(q(1)+q(2)) params(7).*cos(q(1)+q(2)); 0 0];

%Calcolo del vettore dei termini gravitazionali:
g = -[(params(1).*transpose(g0)*Jp1(:,1) + params(5).*transpose(g0)*Jp2(:,1)); (params(1).*transpose(g0)*Jp1(:,2) + params(5).*transpose(g0)*Jp2(:,2))];

end