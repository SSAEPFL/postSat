% %% Parametres %%
% %%%%%%%%%%%%%%%%
% 
repertoire = ' '; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'a.exe'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.dat'; % Nom du fichier d'entree de base MODIFIER SELON VOS BESOINS
fichier = 'output.out';

%tfin = 3.155814954e7;
%nsteps = [0.5,1, 3, 5,9,15,25,27,45,75,135,225,675,46751,140253,233755]; % Nombre d'iterations entier de 10^2 a 10^4  MODIFIER SELON VOS BESOINS
%nsteps = tfin./nsteps;
% tfin =5579.27;
% dt = [0.01,0.1,0.5,1,2, 4, 5,8,10,16,20,40,71,80,142,284,355,568,710,1136,1420,2840,5680];
% nsteps = tfin./dt;
dr = linspace(0.04273,0.04455,100);
 nsimul = length(dr); % Nombre de simulations a faire
% paramstr = 'nsteps'; % Nom du parametre a scanner  MODIFIER SELON VOS BESOINS
% param = vitesse; % Valeurs du parametre a scanner  MODIFIER SELON VOS BESOINS
paramstr = 'dr'; % Nom du parametre a scanner  MODIFIER SELON VOS BESOINS
param = dr; % Valeurs du parametre a scanner  MODIFIER SELON VOS BESOINS

%% Simulations %% 
%%%%%%%%%%%%%%%%%
% Lance une serie de simulations (= executions du code C++)
% Normalement, on ne devrait pas avoir besoin de modifier cette partie

 output = cell(1, nsimul); % Tableau de cellules contenant les noms des fichiers de sortie
 for i = 1:nsimul
     output{i} = [paramstr, '=', num2str(param(i)), '.out'];
     % Execution du programme en lui envoyant la valeur a scanner en argument
      cmd = sprintf('%s%s %s %s=%.15g output=%s', repertoire, executable, input, paramstr, param(i), output{i});
      disp(cmd);
      system(cmd);
 end
 %% Ref loading
 data_ref = load("output_ref.out");
 %% JE veux trouver le dr qui colle au potentiel ponctuel pour le géopotentiel d'ordre 0
delta = zeros(1,nsimul);

for i = 1:length(dr)
   data = load(output{i});
%     figure
     delta(i) = max(abs(sqrt(data(:,2).^2 + data(:,3).^2 + data(:,4).^2) - sqrt(data_ref(:,2).^2+ data_ref(:,3).^2 + data_ref(:,4).^2)));
end
%% Plot
figure 
plot(dr,delta,'rx-')
xlabel('dr [m]','fontsize',18,'interpreter','latex')
ylabel('$\Delta h$ [m]','fontsize',18,'interpreter','latex')
grid on
