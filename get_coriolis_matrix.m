function [C] = get_coriolis_matrix(q, vq, params)
%Funzione per ricavare il vettore dei termini gravitazionali g.

%q: vettore delle coordinate lagrangiane (theta1 theta2)
%vq: vettore delle velocità delle coordinate lagrangiane (vtheta1 vtheta2)
%params: vettore dei parametri caratteristici del robot (m1 a1 l1 I1 m2 a2 l2 I2)

%Simboli di Christoffel:
h = -params(5).*params(2).*params(7).*sin(q(2));
sc = [0 h; -h 0];
sc(:,:,2) = [h h; 0 0];                 %matrice 3D per raccogliere simboli Christoffel

%Ricavo matrice di Coriolis:
C = zeros(2,2);
for i= 1:2
    for j = 1:2
        for k = 1:2
            C(i,j) = C(i,j) + sc(i,j,k).*vq(k);
        end
    end
end
end