% TEST CBF MANIPOLATORE PLANARE 2 LINKS

clc
clear
close all

dt = 0.01;

% Parametri robot: params = (m1 a1 l1 I1 m2 a2 l2 I2)
params = [1, 1, 0.5, 1, 1, 1, 0.5, 1];

% Condizioni iniziali:
q = [-1;2];                        
dq  = [0; 0];
ddq = [0; 0];

% Guadagni regolatore:
kp = 20;                                    %parte proporzionale
kd = 2*sqrt(kp);                            %parte derivativa

% Posizione desiderata:
qB = [1;1];

% Inizializzo errore:
e = qB - q;                              %errore
de = -dq;                                    %derivata dell'errore

% Definizione punti barriera:                                  
x = [1.8; -4; -3.5; 2];                     %q1
y = [2.9; 3; 0.4; 0.3];                     %q2
c = boundary(x,y);
cl = length(c);                             %n punti poligonale

for nn1 = 1:cl-1
        
    p1 = c(nn1);                                                        %pedice punto 1
    p2 = c(nn1+1);                                                      %pedice punto 2

    h(nn1,1) = (y(p2) - y(p1))/(x(p2) - x(p1));                         %coeff angolare 1
    h(nn1,2) = -1;                                                      %coeff angolare 2
    h(nn1,3) = y(p1) - (y(p2) - y(p1))/(x(p2) - x(p1))*x(p1);           %intercetta asse y

end

figure
plot(x(c),y(c))
grid on
title('Barriera')
xlabel('q1')
ylabel('q2')


% Variabili utili per i cicli:
nmax = 1500;
run_loop = 1;                               
i = 0;  

% Posizionamento robot:
while run_loop && i < nmax
        i = i+1;

        % Input di controllo nominale
        tau = kp.*e + kd.*de;                                            %ingresso di controllo
        [B, C, g] = get_dynamics(q,dq,params);                           %vettore contenente [B,C,g]
        n = C*dq + g;
        torque_control = computed_torque_control(dq,tau,B, C, g);        %vettore coppia di controllo

        % Control Barrier Function
        invB = pinv(B);
        H = 2*eye(2);
        f = - 2 * computed_torque_control(dq,tau,B, C, g)';
        A = zeros(cl-1,2);
        b = zeros(cl-1,1);
        gamma1 = 1;
        gamma2 = 2;
        for ii = 1:cl-1
            mi = h(ii,1:2);
            bi = h(ii,3);
            A(ii,1:2) = mi*invB;
            b(ii) = - gamma1*(mi*q + bi) - gamma2*(mi*dq) + (mi*invB)*n;
        end
        u = quadprog(H, f, A, b, ...
        [], [], [], [], [], optimoptions('quadprog','Display','off'));

        % Integration step
        ddq = pinv(B)*(u-n);
        dq  = dq + ddq*dt;
        q   = q + dq*dt;
    
        % Aggiornamento errore
        e = qB - q;
        de = -dq;

        % Salvataggio dati intermedi
        res_robot(i,:) = q;

        % Condizioni fine ciclo
        if norm(e) < 0.01 && norm(de) < 0.001
            run_loop = 0;
        end
end

figure; 
hold on;
plot(res_robot)
plot(repmat(qB,1,length(res_robot))','r--')
legend('q1','q2','Location','northwest')
grid on
title('Grafico posizionamento robot')
xlabel('n iterazioni')
ylabel('Coordinate lagrangiane [rad]')











