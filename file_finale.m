% TESI MANIPOLATORE PLANARE A DUE BRACCI

close all
clear
clc

% Parametri robot: params = (m1 a1 l1 I1 m2 a2 l2 I2)
params = [1, 1, 0.5, 1, 1, 1, 0.5, 1];

% Condizioni di riposo iniziali:
pA = [2; 0];                               %coordinate punto A
qA = inverse_kinematics(pA, 1, params);    %cinematica inversa punto A

% Condizioni inziali robot simulato:
q = [qA(1); qA(2)];                        
dq  = [0; 0];
ddq = [0; 0];

% Condizioni inziali robot fisico:
q_f = [qA(1); qA(2)];
dq_f  = [0; 0];
ddq_f = [0; 0];

% Inizializzo variabili di errore:
% % e = zeros(2,1);                            %errore
de = -dq;                                  %derivata dell'errore

% Frequenza del controllo:
dt = 0.01;                                 %controllo a 100 Hz

% Guadagni regolatore:
kp = 20;                                   %guadagno regolatore parte proporzionale
kd = 2*sqrt(kp);                           %guadagno regolatore parte derivativa

% Variabili utili per il codice:
nmax = 1500;                               %numero max cicli
npos = 10;                                 %numero posizionamenti
c = zeros;                                 %output boundary

% Inizializzo matrice di memorizzazione delle posizioni:
Mm = zeros(npos,2);                        %matrice contenente q

% Plot ostacolo (cerchio):
r = 1;                                   %raggio
c_o = [-1 2];                             %centro
n_o = 1000;                                %precisione disegno
t = linspace(0,2*pi,n_o);
x_o = c_o(1) + r*sin(t);
y_o = c_o(2) + r*cos(t);
line(x_o,y_o)                              
fill(x_o,y_o,'r')
hold on
axis([-5 3 0 3.5])
axis equal
grid on
title('Grafico posizionamento robot')
xlabel('q1')
ylabel('q2')
legend('ostacolo','location','northwest')

% Ciclo:
for pos = 1:npos
    % Inizializzo punto partenda entrambi robot
    qA = q;
    q_f = qA;
    
    % Inizializzo punto arrivo entrambi robot
    pB = random_pos(params);                            %coordinate punto B
    qB = inverse_kinematics(pB, 1, params);             %cinematica inversa punto B

    % Robot simulato:
    e = qB - qA;                                        
    de = -dq;
    run_loop = 1;                                       %condizione cicli
    i = 0;                                              %indice cicli

    while run_loop && i < nmax
        i = i+1;

        % Dinamica diretta:
        tau = kp.*e + kd.*de;                                            %ingresso di controllo
        [B, C, g] = get_dynamics(q,dq,params);                           %vettore contenente [B,C,g]
        n = C*dq + g;
        torque_control = computed_torque_control(dq,tau,B, C, g);        %vettore coppia di controllo
        invB = pinv(B);

        % Step di integrazione:
        ddq = pinv(B)*(torque_control-n);
        dq  = dq + ddq*dt;
        q   = q + dq*dt;
    
        % Aggiornamento errore:
        e = qB - q;
        de = -dq;
        
        % Condizioni fine ciclo:
        if norm(e) < 0.01 && norm(de) < 0.001
            run_loop = 0;
        end
    end

    % Plot posizione end-effector robot simulato:
    hold on
    scatter(q(1),q(2),'.')
    legend('q','location','northwest')
    pause(2)

    % Costruzione matrice con coordinate lagrangiane:
    Mm(pos,1:2) = q;
    x = Mm(:,1);
    y = Mm(:,2);
    
    % Condizione sulla presenza dell'ostacolo:
    dist = sqrt((q(1) - c_o(1))^2 + (q(2) - c_o(2))^2);
    if dist > r
    % Contorno posizioni raggiunte:
        if pos >= 4
            c = boundary(x,y,0.1);                                                  

            % Costruzione bordi:
            cl = length(c);
            for nn1 = 1:cl-1
        
                p1 = c(nn1);                                                        %pedice punto 1
                p2 = c(nn1+1);                                                      %pedice punto 2

                h(nn1,1) = (y(p2) - y(p1))/(x(p2) - x(p1));                         %coeff angolare 1
                h(nn1,2) = -1;                                                      %coeff angolare 2
                h(nn1,3) = y(p1) - (y(p2) - y(p1))/(x(p2) - x(p1))*x(p1);           %intercetta asse y

            end

            % Robot fisico:
            e = qB - q_f;                                        
            de = -dq_f;
            run_loop = 1;                                       %condizione cicli
            i = 0;                                              %indice cicli

            while run_loop && i < nmax
                i = i+1;

                % Dinamica diretta
                tau = kp.*e + kd.*de;                                            %ingresso di controllo
                [B, C, g] = get_dynamics(q_f,dq_f,params);                           %vettore contenente [B,C,g]
                n = C*dq_f + g;
                torque_control = computed_torque_control(dq_f,tau,B, C, g);        %vettore coppia di controllo
                
                % Control Barrier Function
                invB = pinv(B);
                H = 2*eye(2);
                f = - 2 * computed_torque_control(dq_f,tau,B, C, g)';
                A = zeros(cl-1,2);
                b = zeros(cl-1,1);
        
                alpha1 = 1;      %prova con 0.1, 1, 10
                alpha2 = 1;
                gamma1 = alpha2;
        
                for ii = 1:cl-1
                    mi = h(ii,1:2);
                    bi = h(ii,3);
                    A(ii,1:2) = - mi*invB;
                    b(ii) = gamma1*(mi*dq_f) + alpha2*(mi*dq_f + alpha1*(mi*q_f + bi)) - mi*invB*n;
                    if x(c(ii+1)) > x(c(ii))
                        A(ii,1:2) = - A(ii,1:2);
                        b(ii) = - b(ii);
                    end
                end

                u = quadprog(H, f, A, b, ...
                [], [], [], [], [], optimoptions('quadprog','Display','off'));

                % Step di integrazione:
                ddq_f = pinv(B)*(u-n);
                dq_f  = dq_f + ddq_f*dt;
                q_f   = q_f + dq_f*dt;
    
                % Aggiornamento errore:
                e = qB - q_f;
                de = -dq_f;
        
                % Condizioni fine ciclo:
                if norm(e) < 0.01 && norm(de) < 0.001
                    run_loop = 0;
                end
            end
            hold on
            plot(x(c),y(c))
            legend('barriera','location','northwest')
            pause(2)
            scatter(q_f(1),q_f(2),'*')
            legend('qf','location','northwest')
            pause(2)
            hold off
        else
            pos = pos - 1;
        end
    end
end