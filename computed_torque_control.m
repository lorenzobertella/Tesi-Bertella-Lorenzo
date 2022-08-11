function [tau1] = computed_torque_control(vq,tau,B,C,g)
%Funzione che valuta l'ingresso di controllo

%vq: vettore delle velocità delle coordinate lagrangiane (vtheta1 vtheta2)
%tau: termine ingresso di controllo dipendente dall'errore 
%B: matrice d'inerzia
%C: matrice di Coriolis
%g: vettore dei termini gravitazionali

%Ricavo ingresso di controllo:
%tau1 = tau + g;                          %solo controllore PD
tau1 = B*tau+C*vq+g;                      %con feedback linearization
end