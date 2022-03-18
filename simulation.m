%

data = load('output.out');



figure
hold on

plot3(data(:,1),data(:,2),data(:,3),'r-')
plot3(data(:,8),data(:,9),data(:,10),'b-')

