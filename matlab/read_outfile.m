function [t,delta,omega,Pm,Eap,Vmag,theta,E1,Efd,P3,Temperature] = read_outfile(fname,ps,opt)
% usage: [t,delta,omega,Pm,Eap,Vmag,theta,E1,Efd,P3,Temperature] = read_outfile(fname,ps,opt)

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
Temperature = X(:,ix.x.temp)+20; % the temperature value in X is referenced to 20

% Y vars
Vmag  = Y(:,ix.y.Vmag);
theta = Y(:,ix.y.theta);

