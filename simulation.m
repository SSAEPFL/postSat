%
% 
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

figure

hold on

plot3(data(:,2)-data(:,5),data(:,3)-data(:,6),data(:,4)-data(:,7),'k-') % Satellite
plot3(data(:,5)-data(:,5),data(:,6)-data(:,6),data(:,7)-data(:,7),'rx') % Terre

xlabel('x')
ylabel('y')
zlabel('z')
