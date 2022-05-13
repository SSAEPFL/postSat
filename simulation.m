data = load('output_terre.out');
data2 = load('output_inertie.out');
data3 = load('output_frot.out');
data4 = load('output_soleil.out');
data5 = load('output_lune.out');
data6 = load('output_pression.out');
data_total = load('output.out');
% figure
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

r = sqrt(data(:,2).^2 + data(:,3).^2 + data(:,4).^2);
r2 = sqrt(data2(:,2).^2 + data2(:,3).^2 + data2(:,4).^2);
r3 = sqrt(data3(:,2).^2 + data3(:,3).^2 + data3(:,4).^2);
r4 = sqrt(data4(:,2).^2 + data4(:,3).^2 + data4(:,4).^2);
r5 = sqrt(data5(:,2).^2 + data5(:,3).^2 + data5(:,4).^2);
r6 = sqrt(data6(:,2).^2 + data6(:,3).^2 + data6(:,4).^2);
figure
%plot(data(:,1), r, 'k-')
hold on
plot(data2(:,1),(r2-r)./r, 'r.-')
%plot(data(:,1),6.371009e6*ones(size(data(:,1))),'r--')
xlabel('$t$ [s]')
ylabel('$r$ [m]')
legend('Inertial forces', 'Earth radius')
figure
plot(data2(:,1),(r3-r)./r, 'r.-')
hold on
%plot(data3(:,1),r3, 'b--')
%plot(data(:,1),6.371009e6*ones(size(data(:,1))),'r--')
xlabel('$t$ [s]')
ylabel('$r$ [m]')
legend('Drag Force', 'Earth radius')
figure 
%plot(data(:,1), r, 'k-')
hold on 
plot(data2(:,1),(r4-r)./r, 'r.-')
%plot(data(:,1),6.371009e6*ones(size(data(:,1))),'r--')
xlabel('$t$ [s]')
ylabel('$r$ [m]')
legend('Suns gravit. force', 'Earth radius')

figure 
plot(data2(:,1),(r5-r)./r, 'r.-')
hold on 
%plot(data4(:,1),r5, 'r.-')
%plot(data(:,1),6.371009e6*ones(size(data(:,1))),'r--')
xlabel('$t$ [s]')
ylabel('$r$ [m]')
<<<<<<< HEAD
legend('Moons gravit. force', 'Earth radius')

figure 
plot(data2(:,1),(r6-r)./r, 'r.-')
hold on 

%plot(data(:,1),6.371009e6*ones(size(data(:,1))),'r--')
xlabel('$t$ [s]')
ylabel('$r$ [m]')
legend('radiation pressure force', 'Earth radius')


r_total = sqrt(data_total(:,2).^2 + data_total(:,3).^2 + data_total(:,4).^2);
figure
plot(data(:,1), r, 'k-')
hold on
plot(data2(:,1),r_total, 'r.-')
=======
>>>>>>> 43767affbeab3b440be841a92ed0c11e59c90a59
