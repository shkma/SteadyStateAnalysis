%% Note
% Program written by: Shoma Kitayama (shomakit@buffalo.edu)
% Last modified: 30.Aug.2020
% Program last used: MATLAB 9.8.0.1380330 (R2020a) Update 2
% Analysis is conducted to solve w/ws for given values of U0/Uy
% i.e., solution of w/ws given a particular U0/Uy

%% Analyzing

clear, clc

% Structural parameters -----
beta_v = 0.1;     % Inherent damping ratio (see Eq.24 in paper)
alpha  = 0.1;     % Stiffness change (see Eq.3 in paper)
AgAy   = 0.6;     % Acceleration ratio, Ag/Ay (see Eq.29, 30 in paper)

% Analysis parameters -----
delta_U0Uy   = 0.01;   % Step of U0/Uy (will be inscreased incrementally)
Max_U0Uy     = 10;     % Maximum U0/Uy considered in analysis
Num_analysis = round(Max_U0Uy/delta_U0Uy); % Number of ananlysis conducted

% Initialize matrix that are filled with steady-state response analysis results -----
wws_one_solution         = zeros(Num_analysis,2); % Size unknown yet so initialize for maximum considered
wws_two_solution_Larger  = zeros(Num_analysis,2); % Size unknown yet so initialize for maximum considered
wws_two_solution_Smaller = zeros(Num_analysis,2); % Size unknown yet so initialize for maximum considered
i_count=0; i_count_LS=0;  % Initialize counting

% Initialize matrix that are filled with stability analysis results -----
wws_one_solution_stab           = zeros(Num_analysis,2); % Size unknown yet so initialize for maximum considered
wws_two_solution_Larger_stab    = zeros(Num_analysis,2); % Size unknown yet so initialize for maximum considered
wws_two_solution_Smaller_stab   = zeros(Num_analysis,2); % Size unknown yet so initialize for maximum considered
i_count_stab=0; i_count_LS_stab = 0;  % Initialize counting

for i1 = 1 : Num_analysis
  
  % Display the progress of processing data
    fprintf('~ Progress: %s th analysis out of %s analysis is currently being processed \n', num2str(i1), num2str(Num_analysis));
    
    U0Uy = delta_U0Uy * i1;  % U0/Uy.. Increasing with small steps (each step: delta_U0Uy)
    
    wws = sym('wws','positive'); % wws=w/ws (see, for example, Eq.29 in paper)
        
        Theta_SC  = acos(1.0/U0Uy); % Theta for U0/Uy (see Eq.15 in paper)
        
        if (U0Uy > 1.0)    % (see Eq.13, 14 in paper)
            C_SC  = 1.0/pi * ( sin(2.0*Theta_SC) - 2.0 * Theta_SC ) + 1.0;
            S_SC  = 0.0;
            
            % - - -  Stability analysis  - - -
            Csc_Prime  = 0.6366*cos(Theta_SC)*(cos(2.0*Theta_SC)-1.0)/sqrt(1.0-(cos(Theta_SC))^2.0);
            Eq_46_left = (-(wws^2.0)+alpha+(1.0-alpha)*Csc_Prime+(1.0-alpha)*C_SC)*(-(wws^2.0)+alpha+(1.0-alpha)*C_SC)...
                          + (2.0*beta_v*wws-(1.0-alpha)*S_SC)^2.0;
            Computed_wws_Stab = solve(Eq_46_left, wws);  % Solve Eq.46 w.r.t. w/ws
                
                % Judgement of number of solution (0 or 1 or 2)
                if size(Computed_wws_Stab) == [0 0]      % If there is 0 solution..
                    % do nothing (no solution to be used, either U/Uy is too large or too small)
                elseif size(Computed_wws_Stab) == [1 1]  % If there is 1 solution..
                    i_count_stab = i_count_stab+1;
                    wws_one_solution_stab(i_count_stab,:)            = double([Computed_wws_Stab(1)  U0Uy]);
                elseif size(Computed_wws_Stab) == [2 1]  % If there is 2 solutions..
                    i_count_LS_stab = i_count_LS_stab+1;
                    wws_two_solution_Larger_stab(i_count_LS_stab,:)  = double([Computed_wws_Stab(2) U0Uy]);  % Large w/ws
                    wws_two_solution_Smaller_stab(i_count_LS_stab,:) = double([Computed_wws_Stab(1) U0Uy]);  % Small w/ws
                end
            % - - -  Stability analysis end  - - -
                
        else    % (see Eq.13, 14 in paper)
            C_SC  = 1.0;
            S_SC  = 0.0;
        end
        
        % Eq.29 in paper
        Eq_29_left  = U0Uy  - AgAy/sqrt( ( (1-alpha)*C_SC + alpha - wws^2)^2 + ( (1-alpha)*S_SC - 2.0*beta_v*wws )^2 );
        Computed_wws = solve(Eq_29_left, wws);
        
        % Judgement of number of solution (0 or 1 or 2)
        if size(Computed_wws) == [0 0]      % 0 solution
           % do nothing (no solution to be used, either U/Uy is too large or too small)
        elseif size(Computed_wws) == [1 1]  % 1 solution
           i_count = i_count+1;
           wws_one_solution(i_count,:) = double([Computed_wws(1) U0Uy]);
        elseif size(Computed_wws) == [2 1]  % 2 solutions
           i_count_LS = i_count_LS+1;
           wws_two_solution_Larger(i_count_LS,:)  = double([Computed_wws(2) U0Uy]);   % Larger w/ws (among obtained two solutions)
           wws_two_solution_Smaller(i_count_LS,:) = double([Computed_wws(1) U0Uy]);   % Smaller w/ws (among obtained two solutions)
        end

end

%% Outputting (writing on file)

% - - - - - - - - - - - - - - - - - - - -
% Output steady-state analysis  - - - - -
% - - - - - - - - - - - - - - - - - - - -
wws_one_solution         = wws_one_solution(1:i_count,:);            % Cut off zero lines
wws_two_solution_Larger  = wws_two_solution_Larger(1:i_count_LS,:);  % Cut off zero lines
wws_two_solution_Smaller = wws_two_solution_Smaller(1:i_count_LS,:); % Cut off zero lines

Data_SteadyState = [wws_one_solution; wws_two_solution_Larger; flipud(wws_two_solution_Smaller(1:end,:))]; % Store data

dlmwrite(strcat('Data_SteadyState', '.txt'), Data_SteadyState, 'delimiter', ' ', 'precision', 8'); % Output data

% - - - - - - - - - - - - - - - - - - - -
% Output stability loci-path  - - - - - -
% - - - - - - - - - - - - - - - - - - - -
wws_one_solution_stab         = wws_one_solution_stab(1:i_count_stab,:);            % Cut off zero lines
wws_two_solution_Larger_stab  = wws_two_solution_Larger_stab(1:i_count_LS_stab,:);  % Cut off zero lines
wws_two_solution_Smaller_stab = wws_two_solution_Smaller_stab(1:i_count_LS_stab,:); % Cut off zero lines

Data_Stability = [wws_one_solution_stab; wws_two_solution_Larger_stab; flipud(wws_two_solution_Smaller_stab(1:end,:))]; % Store data

dlmwrite(strcat('Data_Stability', '.txt'), Data_Stability, 'delimiter', ' ', 'precision', 8'); % Output data

%% Plotting

         % Plot results of steady-state response analysis -----
         NumSubFigRow_1=1; NumSubFigCol_1=2; NumID_1 = 1;
         subplot(NumSubFigRow_1, NumSubFigCol_1, NumID_1); % Position of each singl figures (4*2)

             plot(Data_SteadyState(:,1), Data_SteadyState(:,2), 'r', 'LineWidth',  0.5);
             hold on
             title(strcat('Output steady-state response analysis results'),'FontName','Times New Roman','Fontsize',9);
             xlabel('\omega/\omega_s','Fontname','Times New Roman','Fontsize',10);
             ylabel('\itU \rm_0 / \itU \rm_y','Fontname','Times New Roman','Fontsize',10);
             set(gca,'Fontname','Times New Roman','Fontsize',10);
             set(gcf, 'Position',  [100, 100, 1000, 400]); % Figure size
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

         % Plot results of stability analysis -----
         NumSubFigRow_2=1; NumSubFigCol_2=2; NumID_2 = 2;
         subplot(NumSubFigRow_2, NumSubFigCol_2, NumID_2); % Position of each singl figures (4*2)

             plot(Data_Stability(:,1), Data_Stability(:,2), 'r', 'LineWidth',  0.5);
             hold on
             title(strcat('Output stability analysis results'),'FontName','Times New Roman','Fontsize',9);
             xlabel('\omega/\omega_s','Fontname','Times New Roman','Fontsize',10);
             ylabel('\itU \rm_0 / \itU \rm_y','Fontname','Times New Roman','Fontsize',10);
             set(gca,'Fontname','Times New Roman','Fontsize',10);
             set(gcf, 'Position',  [100, 100, 1000, 400]); % Figure size
             legend('Stability analysis results')
             xlim([0.0 2.0]);
             ylim([0.0 8.0]);
             xticks([0.0 0.5 1.0 1.5 2.0]);
             xticklabels([0.0 0.5 1.0 1.5 2.0]);
             yticks(     [0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0]);
             yticklabels([0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0]);
             ytickformat('%.1f'); % determine the decimal values for y-axis
             grid on
             hold off
