%

data = load('output.out');



figure
title('Soleil 1 an')
hold on

plot3(data(:,2),data(:,3),data(:,4),'rx') % Soleil 1 ans
plot3(data(:,8),data(:,9),data(:,10),'ko', 'Markersize', 18)
xlabel('x')
ylabel('y')
zlabel('z')

figure
title('Lune 28 jours')
hold on

plot3(data(1:1:4483840,5),data(1:1:4483840,6),data(1:1:4483840,7),'rx') % Lune 1 mois
plot3(data(1:1:4483840,8),data(1:1:4483840,9),data(1:1:4483840,10),'ko', 'Markersize', 18)
xlabel('x')
ylabel('y')
zlabel('z')
