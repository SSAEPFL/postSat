data = load('output.out');
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
 plot3(data(:,11),data(:,12),data(:,13),'r-')
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
figure
plot(data(:,1), r, 'k-')
hold on
plot(data(:,1),6.371009e6*ones(size(data(:,1))),'r--')
xlabel('$t$ [s]')
ylabel('$r$ [m]')
