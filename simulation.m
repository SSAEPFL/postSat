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
title('Lune 28 jours')
hold on

plot3(data(:,5),data(:,6),data(:,7),'rx') % Lune 

xlabel('x')
ylabel('y')
zlabel('z')
