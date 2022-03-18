%

data = load('output.out');



figure
plot3(data(:,2),data(:,3),data(:,4))
hold on
plot3(data(:,5),data(:,6),data(:,7))