function ps = case2_ps
% a 2-bus test case

ps.baseMVA = 100.000000;

ps.bus = [...
%ID type Pd Qd Gs Bs area Vmag Vang zone basekV Vmax Vmin lam P lam Q mu Vx mu Vn locX locY 
 1 2   0 0 0 0 1 1.05 0 230 1 1.05 0.95 0 0 0 0 1 2;
 2 Inf 0 0 0 0 1 1.00 0 230 1 1.05 0.95 0 0 0 0 2 2;
];

ps.branch = [...
%from to R X B rateA rateB rateC tap shift status 
 1 2 0.0 0.2 0.0 40 40 40 1 0 1 28.6897 -15.4187 -27.7847 12.8185 0 0 0.310195 0.29142 0 0 0;
 1 2 0.0 0.2 0.0 40 40 40 1 0 1 28.6897 -15.4187 -27.7847 12.8185 0 0 0.310195 0.29142 0 0 0;
];

ps.gen = [...
%bus Pg Qg Qmax Qmin Vsp mBase status Pmax Pmin mu_Px mu_Pn mu_Qx mu_Qn type cost 
 1 100 0 100 -100 1.05 100 1 0 0 0 0 0 0 2 0 1;
];

ps.shunt = [...
%bus P Q frac_S frac_Z status type value 
% 2 70 70 1 0 1 1 1000 0;
];

ps.mac = [...
%gen r Xd  Xdp Xdpp Xq  Xqp Xqpp D   M Ea  Eap Pm delta omega
 1   0 1.0 0.2   0  0.6  0   0   1.5 3 1.2 0   0  0     0;
];
