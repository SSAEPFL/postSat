%

data = load('output.out');



figure
hold on

plot3(data(:,5),data(:,6),data(:,7),'r-')
plot3(data(:,8),data(:,9),data(:,10),'b-')

