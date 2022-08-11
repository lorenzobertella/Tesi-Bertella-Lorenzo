function [output] = inverse_kinematics(P,pcg,params)
%Cinematica inversa: dato un punto ricavare le coordinate lagrangiane

%P: vettore posizione punto P [xp; yp]
%pcg: posizione nella circonferenza goniometrica della direzione del punto,
%1 parte superiore, 2 parte inferiore
%params: vettore dei parametri caratteristici del robot (m1 a1 l1 I1 m2 a2 l2 I2)

%Definisco variabili utili:
xp = P(1);
yp = P(2);
a1 = params(2);
a2 = params(6);

%Ricavo l'angolo theta2:
c2 = (xp.^2 + yp.^2 - a1.^2 - a2.^2)/(2 .* a1.^2 .* a2.^2);     %cos(theta2)
if pcg == 1                 
    s2 = sqrt(1 - c2.^2);                                       %sin(theta2)      
else 
    s2 = - sqrt(1 - c2.^2);
end

theta2 = atan2(s2,c2);

%Ricavo l'angolo theta1:
xs = a1 + a2 .* c2;                                             %punto fittizio
ys = a2 .* s2;

beta = atan2(ys,xs);                             %angolo rotazione sdr fittizio
gamma = atan2(yp,xp);                            %angolo rotazione sdr punto P
theta1 = gamma - beta;

output = [theta1 ; theta2];
end