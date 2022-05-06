data = load('output.out');
modele = load('modeleRK4.out');
% Un tour en 5579.2552
temps = data(:,1);
x_satellite = data(:,2);
y_satellite = data(:,3);
z_satellite = data(:,4);
N = length(x_satellite);
x_moon = data(:,8);
y_moon = data(:,9);
z_moon = data(:,10);

x_satellite_modele = modele(1:N,2);
y_satellite_modele = modele(1:N,3);
z_satellite_modele = modele(1:N,4);
% On va comparer l'influence d'une force sur notre satellite.
%Distance résiduelle eps = 0.380312591663274 mètres (Norme)
% Sans la lune on retourne au point de départ (à eps près)

ecartx = x_satellite-x_satellite_modele;
ecarty = y_satellite-y_satellite_modele;
ecartz = z_satellite-z_satellite_modele;
ecart = sqrt(ecartx.^2 + ecarty.^2 + ecartz.^2);
% On enlève le résidu
ecart = ecart - 0.380312591663274;
figure 
plot(temps(1:1000:end), ecart(1:1000:end),'r-')
xlabel('temps [s]')
ylabel('ecart')
figure
plot3(x_satellite(1:1000:end),y_satellite(1:1000:end),z_satellite(1:1000:end), 'k-')
hold on 
plot3(x_satellite_modele(1:1000:end),y_satellite_modele(1:1000:end),z_satellite_modele(1:1000:end), 'r--')

