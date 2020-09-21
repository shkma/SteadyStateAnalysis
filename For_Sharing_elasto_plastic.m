%% Note
% Program written by: Shoma Kitayama (shomakit@buffalo.edu)
% Last modified: 30.Aug.2020
% Program last used: MATLAB 9.8.0.1380330 (R2020a) Update 2
% Analysis is conducted to solve w/wBi for given values of U0/Uy
% i.e., solution of w/wBi given a particular U0/Uy

%% Analyzing

clear, clc

% Structural parameters -----
beta_v = 0.1;     % Inherent damping ratio (see Eq.58 in paper)
AgAy   = 1.0;     % Acceleration ratio, Ag/Ay (see Eq.55, 56 in paper)
eta    = 0.2;     % Ratio of preload of SC & yield strength of Bi --- eta=FySC/FyBi, see Eq.49 in paper
pse    = 5.0;     % Ratio of stiffness between SC & Bi --- pse=kSC/kBi, see Eq.50 in paper

% Analysis parameters -----
delta_U0Uy   = 0.01;   % Step of U0/Uy (will be inscreased incrementally)
Max_U0Uy     = 10;     % Maximum U0/Uy considered in analysis
Num_analysis = round(Max_U0Uy/delta_U0Uy); % Number of ananlysis conducted

% Initialize matrix that are filled with steady-state response analysis results -----
wwBi_one_solution         = zeros(Num_analysis,2); % Size unknown yet so initialize for maximum considered
wwBi_two_solution_Larger  = zeros(Num_analysis,2); % Size unknown yet so initialize for maximum considered
wwBi_two_solution_Smaller = zeros(Num_analysis,2); % Size unknown yet so initialize for maximum considered
i_count=0; i_count_LS=0;  % Initialize counting

for i1 = 1 : Num_analysis
  
  % Display the progress of processing data
    fprintf('~ Progress: %s th analysis out of %s analysis is currently being processed \n', num2str(i1), num2str(Num_analysis));
    
    U0Uy = delta_U0Uy * i1;  % U0/Uy.. Increasing with small steps (each step: delta_U0Uy)
    
    wwBi = sym('wwBi','positive'); % wwBi=w/wBi (see, for example, Eq.55 in paper)

        Theta_Bi  = acos(1.0-2.0*(1.0/U0Uy));  % Theta for Bi-linear
        Theta_SC  = acos(eta/pse*(1.0/U0Uy));  % Theta for Self-Centering
        
        % Define hysteresis of Bi-linear
        if (U0Uy > 1.0)
            C_Bi  =  1.0/pi * (Theta_Bi - sin(2.0*Theta_Bi)/2.0);
            S_Bi  = -1.0/pi * (sin(Theta_Bi))^2.0;
        else
            C_Bi  = 1.0;
            S_Bi  = 0.0;
        end

        % Define hysteresis of Self-Centering
        if (U0Uy * (pse/eta) > 1.0)
            C_SC  = 1.0/pi * ( sin(2.0*Theta_SC) - 2.0 * Theta_SC ) + 1.0;
            S_SC  = 0.0;
        else
            C_SC  = 1.0;
            S_SC  = 0.0;
        end
        
        % Eq.55 in paper
        Eq_55_left  = U0Uy  - AgAy/sqrt( ( -(wwBi^2)+(C_Bi+pse*C_SC) )^2 + ( -2.0*beta_v*wwBi + (S_Bi + pse*S_SC) )^2 );
        Computed_wwBi = solve(Eq_55_left, wwBi);
        
        % Judgement of number of solution (0 or 1 or 2)
        if size(Computed_wwBi) == [0 0]      % 0 solution
           % do nothing (no solution to be used, either U/Uy is too large or too small)
        elseif size(Computed_wwBi) == [1 1]  % 1 solution
           i_count = i_count+1;
           wwBi_one_solution(i_count,:) = double([Computed_wwBi(1) U0Uy]);
        elseif size(Computed_wwBi) == [2 1]  % 2 solutions
           i_count_LS = i_count_LS+1;
           wwBi_two_solution_Larger(i_count_LS,:)  = double([Computed_wwBi(2) U0Uy]);   % Larger w/wBi (among obtained two solutions)
           wwBi_two_solution_Smaller(i_count_LS,:) = double([Computed_wwBi(1) U0Uy]);   % Smaller w/wBi (among obtained two solutions)
        end

end

%% Outputting (writing on file)

% - - - - - - - - - - - - - - - - - - - -
% Output steady-state analysis  - - - - -
% - - - - - - - - - - - - - - - - - - - -
wwBi_one_solution         = wwBi_one_solution(1:i_count,:);            % Cut off zero lines
wwBi_two_solution_Larger  = wwBi_two_solution_Larger(1:i_count_LS,:);  % Cut off zero lines
wwBi_two_solution_Smaller = wwBi_two_solution_Smaller(1:i_count_LS,:); % Cut off zero lines

Data_SteadyState = [wwBi_one_solution; wwBi_two_solution_Larger; flipud(wwBi_two_solution_Smaller(1:end,:))]; % Store data

dlmwrite(strcat('Data_SteadyState', '.txt'), Data_SteadyState, 'delimiter', ' ', 'precision', 8'); % Output data

%% Plotting

         % Plot results of steady-state response analysis -----
             plot(Data_SteadyState(:,1), Data_SteadyState(:,2), 'r', 'LineWidth',  0.5);
             hold on
             title(strcat('Output steady-state response analysis results'),'FontName','Times New Roman','Fontsize',9);
             xlabel('\omega/\omega_{Bi}','Fontname','Times New Roman','Fontsize',10);
             ylabel('\itU \rm_0 / \itU \rm_{y,Bi}','Fontname','Times New Roman','Fontsize',10);
             set(gca,'Fontname','Times New Roman','Fontsize',10);
             set(gcf, 'Position',  [100, 100, 450, 400]); % Figure size
             legend('Steady-state response analysis results')
             xlim([0.0 2.0]);
             ylim([0.0 8.0]);
             xticks([0.0 0.5 1.0 1.5 2.0]);
             xticklabels([0.0 0.5 1.0 1.5 2.0]);
             yticks(     [0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0]);
             yticklabels([0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0]);
             ytickformat('%.1f'); % determine the decimal values for y-axis
             grid on
             hold off
