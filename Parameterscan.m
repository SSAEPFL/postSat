% % Ce script Matlab automatise la production de resultats
% % lorsqu'on doit faire une serie de simulations en
% % variant un des parametres d'entree.
% % 
% % Il utilise les arguments du programme (voir ConfigFile.h)
% % pour remplacer la valeur d'un parametre du fichier d'input
% % par la valeur scannee.
% %     
% 
% %% Parametres %%
% %%%%%%%%%%%%%%%%
% 
repertoire = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'a.out'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.dat'; % Nom du fichier d'entree de base MODIFIER SELON VOS BESOINS
fichier = 'output.out';

tfin = 3.155814954e7;
nsteps = [0.5,1, 3, 5,9,15,25,27,45,75,135,225,675,46751,140253,233755]; % Nombre d'iterations entier de 10^2 a 10^4  MODIFIER SELON VOS BESOINS
nsteps = tfin./nsteps;
nsimul = length(nsteps); % Nombre de simulations a faire
paramstr = 'nsteps'; % Nom du parametre a scanner  MODIFIER SELON VOS BESOINS
param = nsteps; % Valeurs du parametre a scanner  MODIFIER SELON VOS BESOINS
% 
%% Simulations %% 
%%%%%%%%%%%%%%%%%
% Lance une serie de simulations (= executions du code C++)
% Normalement, on ne devrait pas avoir besoin de modifier cette partie

 output = cell(1, nsimul); % Tableau de cellules contenant les noms des fichiers de sortie
 for i = 1:nsimul
     output{i} = [paramstr, '=', num2str(param(i)), '.out'];
%     % Execution du programme en lui envoyant la valeur a scanner en argument
%     cmd = sprintf('%s%s %s %s=%.15g output=%s', repertoire, executable, input, paramstr, param(i), output{i});
%     disp(cmd);
%     system(cmd);
 end

% %% Analyse %%
% %%%%%%%%%%%%%
% % Ici, on aimerait faire une etude de convergence: erreur fonction de dt, sur diagramme log-log.
% % A MODIFIER ET COMPLETER SELON VOS BESOINS
% 
 error = zeros(1,nsimul);
    xfinal = zeros(1,nsimul);
    yfinal = zeros(1,nsimul);
    zfinal = zeros(1,nsimul);
% 
for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
    dt(i) = tfin./nsteps(i); % Delta t = t_2-t_1
   %
    x_th = data(1,2)-data(1,5); % TODO: Entrer la vraie solution analytique a tfin
    y_th = data(1,3)-data(1,6);
    z_th = data(1,4)-data(1,7); % TODO: Entrer la vraie solution analytique a tfin
    error(i) = sqrt((data(end,2)-data(end,5)-x_th).^2+(data(end,3)-data(end,6)-y_th).^2+(data(end,4)-data(end,7)-z_th).^2); % erreur sur la position finale
    xfinal(i) = data(end,2);
    yfinal(i) = data(end,3);
    zfinal(i) = data(end,4);
end

 coeff_fitx = polyfit(dt(1:13),xfinal(1:13), 1)
  xFitx = linspace(min(dt), max(dt), 100000);
   coeff_fity = polyfit(dt(1:13),yfinal(1:13), 1)
  xFity = linspace(min(dt), max(dt), 100000);
   coeff_fitz = polyfit(dt(1:13),zfinal(1:13), 1)
  xFitz = linspace(min(dt), max(dt), 100000);
  %yFit = polyval(coeff_fit, xFit); 
% coeff_fit_vit = polyfit(log(nsteps(1:10)),log(error2(1:10)), 1)
%  xFit_vit = linspace(min(log(nsteps(1:10))), max(log(nsteps(1:10))), 1000);
%  yFit_vit = polyval(coeff_fit_vit, xFit_vit);
%  coeff_fit_emec = polyfit(log(nsteps(1:10)),log(error3(1:10)), 1)
%  xFit_emec = linspace(min(log(nsteps(1:10))), max(log(nsteps(1:10))), 1000);
%  yFit_emec = polyval(coeff_fit_emec, xFit_emec); 
% 
for i = 1:nsimul
    data = load(output{i});
    error(i) = sqrt((data(end,2)-coeff_fitx(2)).^2+(data(end,3)-coeff_fity(2)).^2+(data(end,4)-coeff_fitz(2)).^2); % erreur sur la position finale
end
 figure
 loglog(dt, error, 'kx-', 'Linewidth', 1.5)

 hold on 

  xlabel('nsteps','Fontsize', 15)
  ylabel('error of the last position [m]','Fontsize', 15)
 %legend('Mesures experimental','Fit : y(x) = -4.15x+26.31', 'Location', 'NorthWest', 'Fontsize', 15)
grid on
% 




