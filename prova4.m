% TESI MANIPOLATORE PLANARE A DUE BRACCI

close all
clear
clc

% Parametri robot: params = (m1 a1 l1 I1 m2 a2 l2 I2)
params = [1, 1, 0.5, 1, 1, 1, 0.5, 1];

% Punto di partenza:
qA = [0;0];
hold on
%scatter(qA(1),qA(2),500,'.')
%   pause(1)

% Condizioni inziali robot simulato:
q = [qA(1); qA(2)];                        
dq  = [0; 0];
ddq = [0; 0];

% Condizioni inziali robot fisico:
q_f = [qA(1); qA(2)];
dq_f  = [0; 0];
ddq_f = [0; 0];

% Frequenza del controllo:
dt = 0.01;                          %controllo a 100 Hz

% Guadagni regolatore:
kp = 20;                            %guadagno regolatore parte proporzionale
kd = 2*sqrt(kp);                    %guadagno regolatore parte derivativa

% Variabili utili per il codice:
nmax = 1500;                        %numero max cicli (robot)
npos = 1000;                        %numero max posizionamenti
c = zeros;                          %output boundary
q_o = [qA(1); qA(2)];               %variabile che memorizza la posizione
num_p = 0;                          %numero punti in più per allargare zona sicura

% Definizione ostacolo (cerchio):
r = 0.7;                            %raggio
c_o = [1 3];                        %centro
n_o = 1000;                         %precisione disegno
t = linspace(0,2*pi,n_o);
x_o = c_o(1) + r*sin(t);
y_o = c_o(2) + r*cos(t);

% Plot ostacolo:
line(x_o,y_o)                              
fill(x_o,y_o,'r')
hold on
axis([-5 3 0 3.5])
axis equal
grid on
title('Grafico posizionamento robot')
xlabel('q1')
ylabel('q2')

% Lista punti che si vogliono raggiungere:
punti = [-2 1;-3 3;-1 3.5;-2 4;-1 2.5;0 3.75;2 3];

% Plotto prima traiettoria:
%hold on
%plot([qA(1), punti(1, 1)],[qA(2), punti(1, 2)],'r--');

% Ciclo:
for pos = 1:npos
    % variabile di passaggio ostacolo
    z = 0;                                     
    % Inizializzo punto partenda entrambi robot
    qA = q_o;
    q = qA;
    q_f = qA;
    
    % Condizione per non eccedere dimensione vettore punti
    if pos > length(punti)
        break
    end

    % Inizializzo punto arrivo entrambi robot
    qB = punti(pos,:)';
    
    % Robot simulato:
    e = qB - qA;                    %errore                                        
    de = -dq;                       %derivata errore
    run_loop = 1;                   %condizione ciclo
    i = 0;                          %indice cicli
    ni = 1;                         %plot posizionamenti intermedi

    while run_loop && i < nmax
        i = i+1;

        % Dinamica diretta:
        %ingresso di controllo:
        tau = kp.*e + kd.*de;       
        %vettore contenente [B,C,g]:
        [B, C, g] = get_dynamics(q,dq,params);                           
        n = C*dq + g;
        %vettore coppia di controllo:
        torque_control = computed_torque_control(dq,tau,B, C, g);       

        % Step di integrazione:
        ddq = pinv(B)*(torque_control-n);
        dq  = dq + ddq*dt;
        q   = q + dq*dt;
    
        % Aggiornamento errore:
        e = qB - q;
        de = -dq;

        % Condizione sulla presenza dell'ostacolo:
        dist = dist2points(q,c_o);
        if dist > r && z == 0
            q_o = q;
        else
            z = 1;
        end
        
%         % Plot traiettoria:
%         if pos >= 2
%             hold on
%             plot([punti(pos-1,1), qB(1)],[punti(pos-1, 2), qB(2)],'r--');
%         end
% 
%         % Plot posizionamenti intermedi robot simulato:
%         if i == ni * 20
%             hold on
%             scatter(q_o(1),q_o(2),100,'.','black')
%             %pause(1)
%             ni = ni + 1;
%         end

        % Condizioni fine ciclo:
        if norm(e) < 0.01 && norm(de) < 0.001
            run_loop = 0;
        end
    end

    % Plot posizionamento robot simulato:
    hold on
    scatter(q_o(1),q_o(2),500,'.')
    pause(1)
               
    if z == 1
        % Incremento punti per zona sicura:
        k = 0.1;
        Ip(1,1:2) = [q_o(1)+0.001, q_o(2) + k];
        Ip(2,1:2) = [q_o(1) + k, q_o(2)+0.001];
        Ip(3,1:2) = [q_o(1)+0.001, q_o(2) - k];
        Ip(4,1:2) = [q_o(1) - k, q_o(2)+0.001];
        for iii = 1:4
            dist = dist2points(Ip(iii,:),c_o);
            if dist > r
                num_p = num_p+1;
                Mm_p(num_p,1:2) = Ip(iii,:);
            end
        end
        
        % Algoritmo aggiramento ostacolo:
        pB = punti(7,:)';               
        vett = [pB(1)-q_o(1); qB(2) - q_o(2)];
        vettn = [-vett(2); vett(1)];
        u_vettn = vettn/norm(vettn);
        a = 0.1;
        q_o1 = q_o + a*u_vettn;
        dist = sqrt((q_o1(1) - c_o(1))^2 + (q_o1(2) - c_o(2))^2);
        if dist < r
           q_o1 = q_o - a*u_vettn;
        end
        punti(pos+1,:) = q_o1';
        punti(pos+2,:) = pB';
    end

    % Memorizzo posizionamenti:
    Mm(pos,1:2) = q_o;
    if num_p ~= 0
        for ind = pos+1:pos+num_p
            Mm(ind,1:2) = Mm_p(ind-pos,:);
        end
    end
    x = [-2; -1; 0.544; 0.645; 0.762; 0.893; 1.049; 1.249; 2; 2.5; -2; -3];
    y = [1;2.5; 3.545; 3.609; 3.659; 3.695; 3.703; 3.660; 3; 4.5;4;3];

    if pos>= 3

        % Creazione poligonale:
        c = [1,12,11,10,9,8,7,6,5,4,3,2,1];
        
        % Definizione parametri barriere:
        cl = length(c);
        for nn1 = 1:cl-1
            
            % Pedici punto 1 e 2:
            p1 = c(nn1);                                                        
            p2 = c(nn1+1);                                                   
            
            % Costruzione matrice h:
            h(nn1,1) = (y(p2) - y(p1))/(x(p2) - x(p1));                         %coeff angolare 1
            h(nn1,2) = -1;                                                      %coeff angolare 2
            h(nn1,3) = y(p1) - (y(p2) - y(p1))/(x(p2) - x(p1))*x(p1);           %intercetta asse y

        end

        % Robot fisico:
        if z==1
            qB = q_o;
        end
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
                
            % Control Barrier Function:
            %definizione componenti per il quadprog:
            invB = pinv(B);
            H = 2*eye(2);
            f = - 2 * computed_torque_control(dq_f,tau,B, C, g)';
            A = zeros(cl-1,2);
            b = zeros(cl-1,1);
            %parametri dell'equazione di vincolo:
            alpha1 = 1;      %prova con 0.1, 1, 10
            alpha2 = 1;
            gamma1 = alpha2;
            %ciclo per costruire le n equazioni di vincolo:
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
            %nuovo ingresso di controllo:
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
       pause(1)
       scatter(q_f(1),q_f(2),100,'*')
       pause(1)
       hold off

     end
end
 hold on
% plot(x(c),y(c))
% scatter(q_f(1),q_f(2),'*')
% hold off

