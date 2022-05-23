data = load('output.out');
% title('Ascension ')
% plot(data(:,1),data(:,2),'rx') % Soleil 1 ans
% figure
% title('Declinaison ')
% plot(data(:,1),data(:,3),'rx') % Soleil 1 ans
%
%
% figure
% title('Soleil 1 an')
% hold on
%
% plot3(data(:,2),data(:,3),data(:,4),'r-') % Soleil 1 ans
% plot3(data(:,8),data(:,9),data(:,10),'ko', 'Markersize', 18)
% xlabel('x')
% ylabel('y')
% zlabel('z')
%
distance = zeros(1,(length(data)-1)/2);
for i = 1:(length(data)+1)/2
    distance(i) = sqrt((data(i+(length(data)-1)/2,2)-data(1,2))^2+(data(i+(length(data)-1)/2,3)-data(1,3))^2+(data(i+(length(data)-1)/2,4)-data(1,4))^2);
end
[minus index] = min(distance);
r = sqrt(data(:,2).^2 + data(:,3).^2 + data(:,4).^2);
figure
plot(data(:,1), r, 'k-')
hold on
% plot(data2(:,1),(r2-r)./r, 'r.-')
%plot(data(:,1),6.371009e6*ones(size(data(:,1))),'r--')
xlabel('$t$ [s]')
ylabel('$r$ [m]')
legend('Inertial forces', 'Earth radius')

% Make unit sphere
[x,y,z] = sphere;
% Scale to desire radius.
radius = 6.371009e6;
x = x * radius;
y = y * radius;
z = z * radius;

figure 
%plot3(data(:,5),data(:,6),data(:,7),'rx') % Terre
hold on
 %plot3(data(:,11),data(:,12),data(:,13),'r-')
plot3(data(:,2),data(:,3),data(:,4), 'k-') % Satellite
plot3(data(end,2),data(end,3),data(end,4), 'rx') % Satellite
plot3(data(1,2),data(1,3),data(1,4), 'bx') % Satellite
%plot3(data(:,11),data(:,12),data(:,13), 'b-')% SUn
%plot3(data(i,2),data(i,3),data(i,4), 'k-') % Satellite
%plot3(data(i,8),data(i,9),data(i,10), 'bo') %MOon

surf(x,y,z) % terre
hold off
xlabel('x')
ylabel('y')
zlabel('z')

%% 
data_terre = load('output_terre.out');
r_terre = sqrt(data_terre(:,2).^2 + data_terre(:,3).^2 + data_terre(:,4).^2);
figure
plot(data_terre(:,1), r_terre, 'k-')
hold on
% plot(data2(:,1),(r2-r)./r, 'r.-')
% plot(data_terre(:,1),6.371009e6*ones(size(data_terre(:,1))),'r--')
grid on
distance_devi_terre = (max(max(r_terre(1:465))) -max(max(r_terre(71535:end))))/1000 % en KM
distance_oscillation_terre = (max(max(r_terre)) -min(min(r_terre)))/1000 - distance_devi% en KM
xlabel('$t$ [s]')
ylabel('$r$ [m]')
legend('No ext forces', 'Earth radius')
%% 
data_inertie = load('output_inertie.out');
r_inertie = sqrt(data_inertie(:,2).^2 + data_inertie(:,3).^2 + data_inertie(:,4).^2);
figure
plot(data_inertie(:,1), (r_inertie-r_terre)./r_terre, 'k-')
hold on
% plot(data2(:,1),(r2-r)./r, 'r.-')
% plot(data_inertie(:,1),6.371009e6*ones(size(data_inertie(:,1))),'r--')
grid on
distance_devi_inertie = (max(max(r_inertie(1:465))) -max(max(r_inertie(71535:end))))/1000 % en KM
distance_oscillation_inertie = (max(max(r_inertie)) -min(min(r_inertie)))/1000 - distance_devi% en KM
xlabel('$t$ [s]')
ylabel('$r$ [m]')
legend('Inertial Forces', 'Earth radius')
%% 
data_terre = load('output_terre.out');
r_terre = sqrt(data_terre(:,2).^2 + data_terre(:,3).^2 + data_terre(:,4).^2);
figure
plot(data_terre(:,1), r_terre, 'k-')
hold on
% plot(data2(:,1),(r2-r)./r, 'r.-')
% plot(data_terre(:,1),6.371009e6*ones(size(data_terre(:,1))),'r--')
grid on
distance_devi_terre = (max(max(r_terre(1:465))) -max(max(r_terre(71535:end))))/1000 % en KM
distance_oscillation_terre = (max(max(r_terre)) -min(min(r_terre)))/1000 - distance_devi_terre% en KM
xlabel('$t$ [s]')
ylabel('$r$ [m]')
legend('No ext forces', 'Earth radius')
data_sun = load('output_sun.out');
r_sun = sqrt(data_sun(:,2).^2 + data_sun(:,3).^2 + data_sun(:,4).^2);
figure
plot(data_sun(:,1), (r_sun), 'k-')
hold on
% plot(data2(:,1),(r2-r)./r, 'r.-')
% plot(data_sun(:,1),6.371009e6*ones(size(data_sun(:,1))),'r--')
grid on
distance_devi_sun = (max(max(r_sun(1:465))) -max(max(r_sun(71535:end))))/1000 % en KM
distance_oscillation_sun = (max(max(r_sun)) -min(min(r_sun)))/1000 - distance_devi_sun% en KM
xlabel('$t$ [s]')
ylabel('$r$ [m]')
legend('Suns gravit. force', 'Earth radius')
%% 
data_Moon = load('output_Moon.out');
r_Moon = sqrt(data_Moon(:,2).^2 + data_Moon(:,3).^2 + data_Moon(:,4).^2);
figure
plot(data_Moon(:,1), (r_Moon-r_terre)./r_terre, 'k-')
hold on
% plot(data2(:,1),(r2-r)./r, 'r.-')
% plot(data_Moon(:,1),6.371009e6*ones(size(data_Moon(:,1))),'r--')
grid on
distance_devi_moon= (max(max(r_Moon(1:465))) -max(max(r_Moon(71535:end))))/1000 % en KM
distance_oscillation_moon = (max(max(r_Moon)) -min(min(r_Moon)))/1000 - distance_devi% en KM
xlabel('$t$ [s]')
ylabel('$r$ [m]')

legend('Moons gravit. force', 'Earth radius')
%% 
data_pressure = load('output_pressure.out');
r_pressure = sqrt(data_pressure(:,2).^2 + data_pressure(:,3).^2 + data_pressure(:,4).^2);
figure
plot(data_pressure(:,1), (r_pressure-r_terre)./r_terre, 'k-')
hold on
% plot(data2(:,1),(r2-r)./r, 'r.-')
% plot(data_pressure(:,1),6.371009e6*ones(size(data_pressure(:,1))),'r--')
grid on
distance_devi_pressure = (max(max(r_pressure(1:465))) -max(max(r_pressure(71535:end))))/1000 % en KM
distance_oscillation_pressure = (max(max(r_pressure)) -min(min(r_pressure)))/1000 - distance_devi% en KM
xlabel('$t$ [s]')
ylabel('$r$ [m]')
legend('Radiation pressure', 'Earth radius')
%%
data_frot = load('output_frot.out');
r_frot = sqrt(data_frot(:,2).^2 + data_frot(:,3).^2 + data_frot(:,4).^2);
figure
plot(data_frot(:,1), (r_frot-r_terre)./r_terre, 'k-')
hold on
% plot(data2(:,1),(r2-r)./r, 'r.-')
% plot(data_frot(:,1),6.371009e6*ones(size(data_frot(:,1))),'r--')
grid on
distance_devi_frot = (max(max(r_frot(1:465))) -max(max(r_frot(71535:end))))/1000 % en KM
distance_oscillation_frot = (max(max(r_frot)) -min(min(r_frot)))/1000 - distance_devi% en KM
xlabel('$t$ [s]')
ylabel('$r$ [m]')
legend('Drag Force', 'Earth radius')
%% 
data_total = load('output_total.out');
r_total = sqrt(data_total(:,2).^2 + data_total(:,3).^2 + data_total(:,4).^2);
figure
plot(data_total(:,1), (r_total-r_terre)./r_terre, 'k-')
hold on
% plot(data2(:,1),(r2-r)./r, 'r.-')
% plot(data_total(:,1),6.371009e6*ones(size(data_total(:,1))),'r--')
grid on
distance_devi_total = (max(max(r_total(1:5590*2))) -max(max(r_total(858410:end))))/1000 % en KM
distance_oscillation_total = (max(max(r_total(858421:end))) -min(min(r_total(858421:end))))/1000 % en KM
xlabel('$t$ [s]')
ylabel('$r$ [m]')
legend('Total Forces', 'Earth radius')
%%  Temps de simulation -> Nombre de pas = temps simuler , sampling = 1 
second_simuler = [5580,86400, 864000,2592000,15552000, 3.15581e+07]; %Runge-Kutta
elapsed_time = [0.559116,7.65135,77.7391,226.791,1469.37,3529.63]
