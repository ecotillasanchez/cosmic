function [t,delta,omega,Pm,Eap,Vmag,theta,E1,Efd,P3,Temperature,If,It,Pf,Qf,P_inj2,Q_inj2] = read_outfile_attack(fname,ps,opt,outputs)
% usage: [t,delta,omega,Pm,Eap,Vmag,theta,E1,Efd,P3,Temperature,If,It,Pf,Qf,P_inj2,Q_inj2] = read_outfile_attack(fname,ps,opt,outputs)
% custom outfile setup for load oscillating attacks, as described in 
% F. Alanazi, J. Kim and E. Cotilla-Sanchez, "Load Oscillating Attacks 
% of Smart Grids: Vulnerability Analysis," in IEEE Access, vol. 11, 
% pp. 36538-36549, 2023, doi: 10.1109/ACCESS.2023.3266249.

n  = size(ps.bus,1);
ng = size(ps.gen,1);
m  = size(ps.branch,1);
n_sh = size(ps.shunt,1);
ix   = get_indices(n,ng,m,n_sh,opt);

data = readmatrix(fname,'OutputType','double','FileType','text','NumHeaderLines',1);

t = data(:,1);
X = data(:,1+(1:ix.nx));
Y = data(:,1+ix.nx+(1:ix.ny));

% X vars
delta = X(:,ix.x.delta);
omega = X(:,ix.x.omega)*2*pi*ps.frequency;
Pm    = X(:,ix.x.Pm);
Eap   = X(:,ix.x.Eap);
E1    = X(:,ix.x.E1);
Efd   = X(:,ix.x.Efd);
P3   = X(:,ix.x.P3);
Temperature = X(:,ix.x.temp)+20; % the temerature value in X is reference to 20

% Y vars
Vmag  = Y(:,ix.y.Vmag);
theta = Y(:,ix.y.theta);

V_all = Vmag'.*exp(1i.*theta');



If_all = ps.Yf*V_all; % branch status is accounted for in Yf
It_all = ps.Yt*V_all; % branch status is accounted for in Yt
If = If_all;
It = It_all;

% bus injection calculation
I_inj = ps.Ybus*V_all;
S_inj2 = V_all.*conj(I_inj);
S_inj2 = S_inj2';
P_inj2 = real(S_inj2);
Q_inj2 = imag(S_inj2);
% finding_event = find(outputs.event_record(:,2)==11 | outputs.event_record(:,2)==15);
% event_time = outputs.event_record(finding_event,1);
% event_location = outputs.event_record(finding_event,4);
% 
% for i=1:size(finding_event,1)
%     tripping_time=find(round(t,1)==round(event_time(i),1));
%     If(tripping_time(1):end,event_location(i))=0;
% end

F = ps.bus_i(ps.branch(:,1));
%T = ps.bus_i(ps.branch(:,2));
% V_all = V_all';
Sf = V_all(F) .* conj(If);
%St = V_all(T) .* conj(It);
Pf = real(Sf) * ps.baseMVA;
%Pt = real(St) * ps.baseMVA;
Qf = imag(Sf) * ps.baseMVA;
%Qt = imag(St) * ps.baseMVA;

