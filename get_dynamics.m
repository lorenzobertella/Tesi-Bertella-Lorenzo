function [B, C, g] = get_dynamics(q, vq, params)
%Funzione che calcola la matrice d'inerzia B, la matrice di Coriolis e
%il vettore dei termini gravitazionali g

%q: vettore delle coordinate lagrangiane (theta1 theta2)
%vq: vettore delle velocità delle coordinate lagrangiane (vtheta1 vtheta2)
%params: vettore dei parametri caratteristici del robot (m1 a1 l1 I1 m2 a2 l2 I2)

%Ricavo la matrice B:
B = get_inertia_matrix(q,params);

%Ricavo la matrice C:
C = get_coriolis_matrix(q,vq,params);

%Ricavo il vettore g:
g = get_gravity_vector(q,params);

end