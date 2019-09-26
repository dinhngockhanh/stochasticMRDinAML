function PROJECT_growth_rate
%-------------------------------------------Parameters for the algorithm
%   Number of patients
    N_patients          = 6;
%   Number of trajectories per parameter set
    N_traj              = 10;
%   Number of parameter sets
    NO_PARAMETERS       = 1;
%   Weights of all patients
    weight          = [70.3   87.5   103.6  80.7   90.2   79.3];
%   Constants for two negative feedback systems:
    C1                  = 10^-12;
    C2                  = 1;
    k                   = 10^-12;

%   Number of normal cells in bone marrow at equilibrium
    normal_in_BM    = [4.5*10^11 ...
                       4.6*10^11 ...
                       5.0*10^11 ...
                       5.0*10^11 ...
                       4.8*10^11 ...
                       4.6*10^11];
%   Number of normal cells in blood at equilibrium
    normal_in_blood = [1.8*10^10 ...
                       2.0*10^10 ...
                       2.2*10^10 ...
                       2.0*10^10 ...
                       2.0*10^10 ...
                       2.0*10^10];
%----------------------------------Parameters of clones alive by relapse
%   Number of clones in each patients (not counting the normal clone)
    N_cancer            = zeros(1,N_patients);
%   Proliferation rates of all clones of all patients: (columns are
%   patients, rows belong to different clones)
    p_c                 = zeros(max(N_cancer),N_patients);
    a_c                 = zeros(max(N_cancer),N_patients);
    for i_patient=1:N_patients
        if i_patient==1
            patient_ID          = '400220';
            rel_clonal          = [38   22  40]';
        elseif i_patient==2
            patient_ID          = '426980';
            rel_clonal          = [0    3   0   0   9   88]';
        elseif i_patient==3
            patient_ID          = '452198';
            rel_clonal          = [0    0   0   20  80]';
        elseif i_patient==4
            patient_ID          = '573988';
            rel_clonal          = [15   39  46]';
        elseif i_patient==5
            patient_ID          = '758168';
            rel_clonal          = [39   18  35  8]';
        elseif i_patient==6
            patient_ID          = '804168';
            rel_clonal          = [36   45  19]';
        elseif i_patient==7
            patient_ID          = 'Shlush2-1';
            rel_clonal          = [20.26    3.92    34.80   2.19    8.83    30.00];
        elseif i_patient==8
            patient_ID          = 'Shlush2-2';
            rel_clonal          = [29.09    3.92    25.97   2.19    8.83    30.00];
        elseif i_patient==9
            patient_ID          = 'Shlush3-1';
            rel_clonal          = [3.30     26.54   0.93    4.23    65.00];
        elseif i_patient==10
            patient_ID          = 'Shlush3-2';
            rel_clonal          = [7.53     22.31   0.93    4.23    65.00];
        elseif i_patient==11
            patient_ID          = 'Shlush6';
            rel_clonal          = [0.00     33.51   0.00    9.49    0.00    57.00];
        elseif i_patient==12
            patient_ID          = 'Shlush8';
            rel_clonal          = [82.26    12.74   5.00];
        elseif i_patient==13
            patient_ID          = 'Shlush15';
            rel_clonal          = [0.00     47.47   3.52    10.01   39.00];
        end
        N_clones        = length(rel_clonal)-1;
        filename        = ['parameters_' patient_ID '_centered.txt'];
        if (exist(filename,'file')~=2)
            continue;
        end
        fileID          = fopen(filename,'r');
        formatSpec      = '%f';
        A               = fscanf(fileID,formatSpec);
        fclose(fileID);
        A               = reshape(A,2*N_clones,[])';
        para_a_l        = A(1:N_clones);
        para_p_l        = A(N_clones+1:2*N_clones);
        n_alive_clone   = 0;
        for i=1:length(rel_clonal)-1
            if rel_clonal(i)>0
                n_alive_clone   = n_alive_clone+1;
                p_c(n_alive_clone,i_patient)    = para_p_l(i);
                a_c(n_alive_clone,i_patient)    = para_a_l(i);
            end
        end
        N_cancer(i_patient) = n_alive_clone;
    end
%   Plot background
    figure(1)
    clf
%------------------------------------------------------Find growth rates
    lambda_vec  = zeros(max(N_cancer),N_patients);
    for i=1:N_patients
        % fprintf('Patient %d:\n',i);
        c_1_bar = 2*(10^9)*weight(i);
        s_bar   = 1/(1+k*normal_in_blood(i));
        d_bar   = C1*(normal_in_BM(i)-C2*c_1_bar);
        for j=1:N_cancer(i)
            lambda = 30*(log10(exp(1)))*((2*a_c(j,i)*s_bar-1)*p_c(j,i)-d_bar);
            lambda_vec(j,i) = lambda;
        end
    end
%---------------------------------------------------Find time to relapse
    T_relapse_mean_vec  = zeros(N_patients,1);
    N_count_mean        = 0;
    for patient=1:N_patients
        if patient==1
            patient_ID      = '400220';
            L_induction     = 7;
            L_consolidation = 11;
            L_relapse       = 252;
            N_clones        = 2;
        elseif patient==2
            patient_ID      = '426980';
            L_induction     = 7;
            L_consolidation = 8;
            L_relapse       = 805;
            N_clones        = 5;
        elseif patient==3
            patient_ID      = '452198';
            L_induction     = 7;
            L_consolidation = 11;
            L_relapse       = 505;
            N_clones        = 4;
        elseif patient==4
            patient_ID      = '573988';
            L_induction     = 7;
            L_consolidation = 5;
            L_relapse       = 365;
            N_clones        = 2;
        elseif patient==5
            patient_ID      = '758168';
            L_induction     = 21;
            L_consolidation = 90;
            L_relapse       = 961;
            N_clones        = 3;
        elseif patient==6
            patient_ID      = '804168';
            L_induction     = 7;
            L_consolidation = 11;
            L_relapse       = 235;
            N_clones        = 2;
        elseif patient==7
            patient_ID      = 'Shlush15';
            L_induction     = 7;
            L_consolidation = 30;
            L_relapse       = 300;
            N_clones        = 4;
        end
%       Time of remission
        L_remission         = L_induction+30+L_consolidation;
        N_count_T_relapse   = 0;
        T_relapse_vec       = zeros(1,1);
        for kk=1:3
            for ii=1:N_traj
%               Download the trajectory
                filename    = [patient_ID '/stochastic_' patient_ID '_' num2str(kk) '_' num2str(ii) '.txt'];
                if (exist(filename,'file')~=2)
                    continue;
                end
                fprintf('%d---%d---%d\n',ii,kk,patient);
                fileID      = fopen(filename,'r');
                formatSpec  = '%f';
                A           = fscanf(fileID,formatSpec);
                fclose(fileID);
                A           = reshape(A,1+2+2*N_clones,[])';
%               Output the time and trajectory arrays
                T_entire    = A(:,1);
                if T_entire(end)<L_relapse
                    continue;
                end
                P_entire    = A(:,2:end);
%               Find time to relapse
                found_T_relapse = 0;
                index           = 0;
                while (found_T_relapse==0)
                    index   = index+1;
%                   Find %PB blast
                    PB_blast    = 0;
                    for i_clone=1:N_clones
                        PB_blast    = PB_blast+P_entire(index,2*i_clone+1);
                    end
                    Sum     = PB_blast+P_entire(index,1);
                    PBblast = 100*PB_blast/Sum;
%                   Update if "time to relapse" is found...
                    if ((T_entire(index)>L_remission)&&(PBblast>5))
                        found_T_relapse = 1;
                        T_relapse       = T_entire(index)-L_remission;
                        N_count_T_relapse       = N_count_T_relapse+1;
                        T_relapse_vec(N_count_T_relapse)    = T_relapse;
%                   In case the simulation never reaches relapse...
                    elseif (index==length(T_entire))
                        found_T_relapse = 1;
                    end
                end
            end
        end
        T_relapse_mean      = sum(T_relapse_vec(1:N_count_T_relapse))/N_count_T_relapse;
        T_relapse_mean_vec(patient) = T_relapse_mean;
    end
%-------------------------------------------------------------------Plot
    for i=1:N_patients
        if i==1
            color   = 'r';
        elseif i==2
            color   = 'b';
        elseif i==3
            color   = 'g';
        elseif i==4
            color   = 'm';
        elseif i==5
            color   = 'k';
        elseif i==6
            color   = [0.85 0.85 0.85];
        elseif i==7
            color   = [0.35 0.75 0.95];
        end
        p(i) = scatter(lambda_vec(1,i),T_relapse_mean_vec(i),'o','filled','DisplayName',['Patient ' num2str(i)]); hold on;
        p(i).SizeData           = 200;
        p(i).MarkerFaceColor    = color;
        p(i).LineWidth          = 2;
        p(i).MarkerEdgeColor    = 'k';
        for j=2:N_cancer(i)
            p1 = scatter(lambda_vec(j,i),T_relapse_mean_vec(i),'o','filled'); hold on;
            p1.SizeData         = 200;
            p1.MarkerFaceColor  = color;
            p1.LineWidth        = 2;
            p1.MarkerEdgeColor  = 'k';
        end
    end
    plot([0.3 0.3],[0 800],'--r','LineWidth',2);hold on
    plot([2.0 2.0],[0 800],'--r','LineWidth',2);hold on
%----------------------------------------------------------Find Hill fit
    global vec_growth_rate_Hill vec_T_relapse_Hill
    ii  = 0;
    for i=1:N_patients
        for j=1:N_cancer(i)
            ii  = ii+1;
            vec_growth_rate_Hill(ii)    = lambda_vec(j,i);
            vec_T_relapse_Hill(ii)      = T_relapse_mean_vec(i);
        end
    end

    % [growth_rate_values,T_relapse_values]=Hill_fit;
    % p(N_patients+1) = plot(growth_rate_values,T_relapse_values,'-k','DisplayName','Fitted Hill function'); hold on
    % p(N_patients+1).LineWidth = 2;

    set(gca,'fontsize',20);
    l = xlabel('Growth rate');
    l.FontSize  = 25;
    l = ylabel('Time to relapse (days)');
    l.FontSize  = 25;

    [l,objh]    = legend(p(1:N_patients),'Location','northeast');
    l.FontSize  = 20;
    ch = findobj(objh,'type','patch');
    set(ch,'Markersize',10);

    figure_name = ['plot_growth_rate_vs_time_to_relapse.png'];
    saveas(gcf,figure_name);
end
%=======================================================================
%=======================================================================
%=======================================================================
function [growth_rate_values,T_relapse_values] = Hill_fit
    global vec_growth_rate_Hill vec_T_relapse_Hill

    fun         = @(parameter_set,growth_rate) Hill_func(parameter_set,growth_rate);
    [para,resnorm]= lsqcurvefit(fun,[1600,1,-5,5],vec_growth_rate_Hill,vec_T_relapse_Hill);
    growth_rate_values  = [0:0.01:2];

    A           = para(1);
    B           = para(2);
    C           = para(3);
    n           = para(4);
    T_relapse_values= A./(B+(abs(log10(growth_rate_values)-C)).^n);
    fprintf('A=%f\nB=%f\nC=%f\nn=%f\n',para(1),para(2),para(3),para(4));
end
%=======================================================================
%=======================================================================
%=======================================================================
function T_relapse_output = Hill_func(parameter_set,growth_rate)
    A   = parameter_set(1);
    B   = parameter_set(2);
    C   = parameter_set(3);
    n   = parameter_set(4);
    T_relapse_output    = A./(B+(abs(log10(growth_rate)-C)).^n);
    T_relapse_output(isnan(T_relapse_output))=0.03;
end
