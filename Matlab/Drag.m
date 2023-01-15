%% Data loading
data = load('output.out');
%% Data loading w_D
data_w_D = load('output_w_D.out');
%% Plot parameters
lw = 1;
fs = 15;
fsl = 13;
%% Plot with drag
r = sqrt(data(:,2).^2 + data(:,3).^2 + data(:,4).^2);
figure
plot(data(:,1), r, 'k-','linewidth',lw)
hold on
plot(data(:,1),6.371009e6*ones(size(data(:,1))),'r--','linewidth',lw)
xlabel('$t$ [s]','interpreter','latex','fontsize',fs)
ylabel('$r$ [m]','interpreter','latex','fontsize',fs)
legend('Satellite with drag','Earth','interpreter','latex','fontsize',fsl)
%% Without drag
r_w_D = sqrt(data_w_D(:,2).^2 + data_w_D(:,3).^2 + data_w_D(:,4).^2);
figure
plot(data_w_D(:,1), r_w_D, 'k-','linewidth',lw)
hold on
plot(data_w_D(:,1),6.371009e6*ones(size(data_w_D(:,1))),'r--','linewidth',lw)
xlabel('$t$ [s]','interpreter','latex','fontsize',fs)
ylabel('$r$ [m]','interpreter','latex','fontsize',fs)
legend('Satellite without drag','Earth','interpreter','latex','fontsize',fsl)
%% Altitude
figure
plot(data_w_D(:,1), r_w_D - 6.371009e6*ones(size(data_w_D(:,1))), 'b-','linewidth',lw)
hold on
plot(data(:,1), r-6.371009e6*ones(size(data(:,1))), 'r-','linewidth',lw)
%% Differences
delta_r = r_w_D - r;
figure
plot(data_w_D(:,1), delta_r, 'k-','linewidth',lw)
xlabel('$t$ [s]','interpreter','latex','fontsize',fs)
ylabel('$\Delta r$ [m]','interpreter','latex','fontsize',fs)
legend('Difference of altitude without drag/with drag','interpreter','latex','fontsize',fsl)
