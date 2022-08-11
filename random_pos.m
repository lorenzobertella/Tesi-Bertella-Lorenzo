function [P] = random_pos(params)
%Funzione che genera punti casuali nello spazio di lavoro

%Definizione variabili di interesse:
R = params(2) + params(6);

%Calcolo punto casuale in coordinate polari:
r = R .* rand();
beta = rand() .* 2 .* pi;

%Conversione coordinate polari in cartesiane:
P = [r .* cos(beta); r .* sin(beta)];

end