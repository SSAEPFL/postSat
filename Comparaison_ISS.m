data_ISS = load('ISS.txt');
data = load('output.out');
altitude_ISS = 1.e3*sqrt( data_ISS(:,2).^2 + data_ISS(:,3).^2 + data_ISS(:,4).^2);
altitude = sqrt(data(:,2).^2 + data(:,3).^2 + data(:,4).^2);
temps = 0:360;
temps = temps * 4*60;
figure 
plot(temps,altitude_ISS(1:361),'linewidth',1)
hold on
plot(data(:,1),altitude,'linewidth',1)
xlabel('$t$ [s]', 'interpreter','latex','fontsize',18)
ylabel('$h$[m]', 'interpreter','latex','fontsize',18)
legend('ISS','Simulation $k = 4$','interpreter','latex','fontsize',15)
%plot(data_ISS(1:360,





diff = max(abs(altitude_ISS(1:361) - altitude))
% diff with ponctual 2.7696e+04
% diff with k = 2.7696e+04 check! 1.2114e+05 1.1388e+05 
% 7.8562e+04

%% ISS plot

%% Simulation plot
data = load('output.out');
