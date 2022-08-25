function [dist] = dist2points(P1,P2)
% Distanza tra due punti

dist = sqrt((P1(1) - P2(1))^2 + (P1(2) - P2(2))^2);
end