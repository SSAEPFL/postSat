%% Data loading
data = load('output.out');
data_w_D = load('output_w_D.out');
%% Plot with drag
r = sqrt(data(:,2).^2 + data(:,3).^2 + data(:,4).^2);
figure
plot(data(:,1), r, 'k-')
hold on
plot(data(:,1),6.371009e6*ones(size(data(:,1))),'r--')
xlabel('$t$ [s]')
ylabel('$r$ [m]')

%% Without drag
r_w_D = sqrt(data_w_D(:,2).^2 + data_w_D(:,3).^2 + data_w_D(:,4).^2);
figure
plot(data_w_D(:,1), r_w_D, 'k-')
hold on
plot(data_w_D(:,1),6.371009e6*ones(size(data_w_D(:,1))),'r--')
xlabel('$t$ [s]')
ylabel('$r$ [m]')
%% Differences
delta_r = r_w_D - r;
figure
plot(data_w_D(:,1), delta_r, 'k-')