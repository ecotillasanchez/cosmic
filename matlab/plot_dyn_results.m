function plot_dyn_results(t,xy_out,ix)
% usage: plot_dyn_results(t,xy_out,ix)

x_out = xy_out(1:ix.nx,:);
y_out = xy_out((1:ix.ny)+ix.nx,:);

deltas = x_out(ix.x.delta,:);
omegas = x_out(ix.x.omega,:);

figure(1);
plot(t,deltas);
xlabel('time (sec)');
ylabel('Machine angles (rad)');


figure(2);
plot(t,omegas*60);
xlabel('time (sec)');
ylabel('Machine speeds (Hz)');
