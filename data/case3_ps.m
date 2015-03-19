function ps = case3_ps

ps.casename = 'case3';
ps.baseMVA = 100.000000;

ps.bus = [...
%ID type Pd Qd Gs Bs area Vmag Vang basekV zone Vmax Vmin lam_P lam_Q mu_Vx mu_Vn locX locY 
 1 3 0 0 0 0 1 1.05 0 16.5 1 1.1 0.9 0 0 0 0 0 0;
 2 2 0 0 0 0 1 1 0 16.5 1 1.1 0.9 0 0 0 0 0 0;
 3 1 0 0 0 0 1 1 0 16.5 1 1.1 0.9 0 0 0 0 0 0;
];

ps.branch = [...
%from to R X B rateA rateB rateC tap shift status 
 1 3 0 0.05 0 300 300 300 1 0 1 0 0 0 0 0 0 0 0 0 0 0;
 2 3 0 0.03999 0 300 300 300 1 0 1 0 0 0 0 0 0 0 0 0 0 0;
 1 2 0 0.03999 0 300 300 300 1 0 1 0 0 0 0 0 0 0 0 0 0 0;
];

ps.gen = [...           
%bus Pg Qg Qmax Qmin Vsp mBase status Pmax Pmin mu_Px mu_Pn mu_Qx mu_Qn type cost part_fact RRU RRD id 
 1 71.6 10 1000 -1000 1.04 100 1 1000 10 0 0 0 0 3 5 1 0 0 1;
];

ps.shunt = [...          
%bus P Q frac_S frac_Z status type value 
 2 200 50 1 0 1 1 1000 0;
 3 200 50 1 0 1 1 1000 0;
];

ps.mac = [...
%bus r Xd Xdp Xdpp Xq Xqp Xqpp D M Ea Eap Pm delta_m omega Td0 Td0p 
 1 0 0.146 0.0608 0 0.0969 0.0969 0 23.64 47.28 1.2 0 0 0 0 0 8.96;
];

 
ps.exc=[...
%gen type Ka Ta Tb Ke Te Urmax Urmin Vref Efd E1 
 1 1 0 1 10 100 0.1 5 -5;
];      

ps.gov=[...
%gen type R Tt LCmax LCmin Pmax Pmin Pref 
 1 1 0.05 2 100 -100 2.475 0;
 ];
