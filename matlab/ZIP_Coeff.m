function [ZIP_Coff_P_Q_Final] = ZIP_Coeff()
% Active and Reactive ZIP model per customer (Summer) 

LastName = {'Residential';'';'';'';'';'';'Commercial';'';'';'';'';'';'Industrial'};
Zp = [1.5;1.57;1.56;1.31;0.96;1.18;0.27;0.69;0.77;0.55;0.4;0.76;1.21];
Ip = [-2.31;-2.48;-2.49;-1.94;-1.17;-1.64;-0.33;0.04;-0.84;0.24;-0.4;-0.52;-1.61];
Pp = [1.81;1.91;1.93;1.63;1.21;1.47;1.06;0.27;1.07;0.21;1.01;0.76;1.41];
Zq = [7.41;9.28;10.1;9.2;6.28;8.29;5.48;1.82;8.09;0.55;4.43;6.92;4.35];
Iq = [-11.97;-15.29;-16.75;-15.27;-10.16;-13.67;-9.7;-2.24;-13.65;-0.09;-7.98;-11.75;-7.08];
Pq = [5.55;7.01;7.65;7.07;4.88;6.38;5.22;1.43;6.56;0.54;4.56;5.83;3.72];
ZIP_Coff_P_Q = table(LastName,Zp,Ip,Pp,Zq,Iq,Pq);

% Calculating the ZIP values for the residential customer
mean_Zp_Resid = mean(ZIP_Coff_P_Q.Zp(1:6));

mean_Ip_Resid = mean(ZIP_Coff_P_Q.Ip(1:6));

mean_Pp_Resid = mean(ZIP_Coff_P_Q.Pp(1:6));

mean_Zq_Resid = mean(ZIP_Coff_P_Q.Zq(1:6));

mean_Iq_Resid = mean(ZIP_Coff_P_Q.Iq(1:6));

mean_Pq_Resid = mean(ZIP_Coff_P_Q.Pq(1:6));

% Calculating the ZIP values for the Commercial customer
mean_Zp_Comm = mean(ZIP_Coff_P_Q.Zp(7:12));

mean_Ip_Comm = mean(ZIP_Coff_P_Q.Ip(7:12));

mean_Pp_Comm = mean(ZIP_Coff_P_Q.Pp(7:12));

mean_Zq_Comm = mean(ZIP_Coff_P_Q.Zq(7:12));

mean_Iq_Comm = mean(ZIP_Coff_P_Q.Iq(7:12));

mean_Pq_Comm = mean(ZIP_Coff_P_Q.Pq(7:12));

% The Final ZIP result

LastName2 = {'Residential';'Commercial';'Industrial'};

Resid = [mean_Zp_Resid,mean_Ip_Resid,mean_Pp_Resid,mean_Zq_Resid,mean_Iq_Resid,mean_Pq_Resid];
comm  = [mean_Zp_Comm,mean_Ip_Comm,mean_Pp_Comm,mean_Zq_Comm,mean_Iq_Comm,mean_Pq_Comm];
Inds  = [1.21,-1.61,1.41,4.35,-7.08,3.72];

Zp_2 = [mean_Zp_Resid;mean_Zp_Comm;1.21];
Ip_2 = [mean_Ip_Resid;mean_Ip_Comm;-1.61];
Pp_2 = [mean_Pp_Resid;mean_Pp_Comm;1.41];

Zq_2 = [mean_Zq_Resid;mean_Zq_Comm;4.35];
Iq_2 = [mean_Iq_Resid;mean_Iq_Comm;-7.08];
Pq_2 = [mean_Pq_Resid;mean_Pq_Comm;3.72];
ZIP_Coff_P_Q_Final = table(LastName2,Zp_2,Ip_2,Pp_2,Zq_2,Iq_2,Pq_2);
end