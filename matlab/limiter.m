function [dP22_dP2,P22] = limiter(P2,Pmax,Pmin)
% This function simulates a limiter with Pmin and Pmax as the upper
% and lower limits. The function can be used as a rate limiter. In
% that case, the input should be the derivative of some variable.
% The function has a linear region in which the output is equal to
% the input. There are two small sections which are interpolated.
% They are in regions [Pmin Pmin+0.05*(Pmax-Pmin)] and [Pmin+0.95*(Pmax-Pmin) Pmax]
% P2  : input value to the limiter block
% P22 : output value of the limiter block

Pmax(isinf(Pmax)) = 1e15;        % Put a large number to calculate the upper limit when we have Inf's

P22      = P2;
dP22_dP2 = 1;
a = [-400.0000  %a3;a2;a1;a0     These coefficients are calculated for a limiter with low and high value of 0 and 1
    40.0000                    % a is for the range [0 0.05]
    0
    0];
b = 1.0e+03 * [-0.4;   %b3;b2;b1;b0   b is for the range [0.95 1]
    1.16;
    -1.12;
    0.361];

% DETERMINING THE OUTPUT BASED ON THE INPUT VALUE
set1 = P2 <= Pmin;
P22 (set1) = Pmin (set1);
dP22_dP2 (set1) = 0;

set1 = P2 >= Pmax;
P22 (set1) = Pmax (set1);
dP22_dP2 (set1) = 0;

set1 = P2 >= Pmin + 0.0001*(Pmax - Pmin) & P2 <= Pmin + 0.9999*(Pmax - Pmin);
P22 (set1) = P2(set1);
dP22_dP2 (set1) = 1;

set1 = P2 > Pmin & P2 < Pmin + 0.0001*(Pmax - Pmin);                     % This is the region [Pmin Pmin+0.05*(Pmax-Pmin)], first
x (set1)   = (P2 (set1) - Pmin(set1))./(Pmax(set1) - Pmin(set1));  % the input value is converted to a value for [0 1] limiter
P22 (set1) = a(1)*x(set1).^3+a(2)*x(set1).^2+a(3)*x(set1)+a(4);    % and then the output of this limiter is converted to the output
P22 (set1) = P22 (set1) .*(Pmax(set1) - Pmin(set1)) + Pmin(set1);  % for actual limiter with limits [Pmin Pmax]
dP22_dP2 (set1) = 3*a(1)*x(set1).^2+2*a(2)*x(set1)+a(3);


set1 = P2 < Pmax & P2 > Pmin + 0.9999*(Pmax - Pmin);                     % The region [Pmin+0.95*(Pmax-Pmin) Pmax]
x (set1)   = (P2(set1) - Pmin(set1))./(Pmax(set1) - Pmin(set1));
P22 (set1) = b(1)*x(set1).^3+b(2)*x(set1).^2+b(3)*x(set1)+b(4);
P22 (set1) = P22(set1) .*(Pmax(set1) - Pmin(set1)) + Pmin(set1);
dP22_dP2 (set1)  = 3*b(1)*x(set1).^2+2*b(2)*x(set1)+b(3);


dP22_dP2 = dP22_dP2.';
