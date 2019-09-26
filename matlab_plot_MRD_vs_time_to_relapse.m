function matlab_plotting_MRD_to_relapse_all()
    global MRD_mean_vec T_relapse_mean_vec
%=======================================================================
%   Time of wait after remission to calculate MRD
    L_MRD_wait          = 4*30;
%   Number of trajectories per parameter set
    N_traj              = 1000;
%   Number of patients
    N_patients          = 13;
%   Number of patients for fitting the Hill function
    N_patients_fit      = 6;
%   Number of parameter sets
    NO_PARAMETERS       = 1;
%   Index of patient to exclude for a second Hill function
    i_exclude           = 5;
%=======================================================================
    figure(1);
    clf
    MRD_mean_vec        = zeros(N_patients_fit*3,1);
    T_relapse_mean_vec  = zeros(N_patients_fit*3,1);
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
            patient_ID      = 'Shlush2-1';
            L_induction     = 7;
            L_consolidation = 30;
            L_relapse       = 480;
            N_clones        = 5;
        elseif patient==8
            patient_ID      = 'Shlush2-2';
            L_induction     = 7;
            L_consolidation = 30;
            L_relapse       = 480;
            N_clones        = 5;
        elseif patient==9
            patient_ID      = 'Shlush3-1';
            L_induction     = 7;
            L_consolidation = 10;
            L_relapse       = 180;
            N_clones        = 4;
        elseif patient==10
            patient_ID      = 'Shlush3-2';
            L_induction     = 7;
            L_consolidation = 10;
            L_relapse       = 180;
            N_clones        = 4;
        elseif patient==11
            patient_ID      = 'Shlush6';
            L_induction     = 7;
            L_consolidation = 30;
            L_relapse       = 480;
            N_clones        = 5;
        elseif patient==12
            patient_ID      = 'Shlush8';
            L_induction     = 17;
            L_consolidation = 30;
            L_relapse       = 240;
            N_clones        = 2;
        elseif patient==13
            patient_ID      = 'Shlush15';
            L_induction     = 7;
            L_consolidation = 30;
            L_relapse       = 300;
            N_clones        = 4;
        end
%       Time of remission
        L_remission         = L_induction+30+L_consolidation;
%       Time of MRD calculation
        L_MRD               = L_remission+L_MRD_wait;
%-------
        for kk=1:3
            if kk==1
                direction   = 'left';
            elseif kk==2
                direction   = 'centered';
            elseif kk==3
                direction   = 'right';
            end
            MRD_vec         = zeros(NO_PARAMETERS*N_traj,1);
            N_count_MRD     = 0;
            T_relapse_vec   = zeros(NO_PARAMETERS*N_traj,1);
            N_count_T_relapse= 0;
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
%---------------Find MRD
                index       = find(T_entire>L_MRD,1);
                P_MRD       = P_entire(index,:);
                MRD         = 0;
                for i_clone=1:N_clones
                    MRD     = MRD+P_MRD(2*i_clone+1);
                end
                Sum         = MRD+P_MRD(1);
                MRD         = 100*MRD/Sum;
                N_count_MRD = N_count_MRD+1;
                MRD_vec(N_count_MRD)    = MRD;
%---------------Find time to relapse
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
                        T_relapse   = T_entire(index)-L_remission;
                        N_count_T_relapse   = N_count_T_relapse+1;
                        T_relapse_vec(N_count_T_relapse)    = T_relapse;
%                   In case the simulation never reaches relapse...
                    elseif (index==length(T_entire))
                        found_T_relapse = 1;
                    end
                end
            end
            MRD_mean            = sum(MRD_vec(1:N_count_MRD))/N_count_MRD;
            T_relapse_mean      = sum(T_relapse_vec(1:N_count_T_relapse))/N_count_T_relapse;
            if patient==1
                color   = 'r';
            elseif patient==2
                color   = 'b';
            elseif patient==3
                color   = 'g';
            elseif patient==4
                color   = 'm';
            elseif patient==5
                color   = 'k';
            elseif patient==6
                color   = [0.85 0.85 0.85];
            elseif patient==7
                color   = [0.35 0.75 0.95];
            elseif patient==8
                color   = [0.2 0.25 0.95];
            elseif patient==9
                color   = [0.7 0.01 0.3];
            elseif patient==10
                color   = [0.01 0.6 0.01];
            elseif patient==11
                color   = [0.4 0.4 0.4];
            elseif patient==12
                color   = [0.5 0.2 0.9];
            elseif patient==13
                color   = [0.1 0.2 0.7];
            end
            if patient<=N_patients_fit
                text_legend_name    = ['Patient ' patient_ID ' (Ding et al.)'];
            else
                Shlush_ID = extractAfter(patient_ID,"Shlush");
                text_legend_name    = ['Patient ' Shlush_ID ' (Shlush et al.)'];
            end
            if kk==1
                p1 = scatter(MRD_mean,T_relapse_mean,'<','filled'); hold on;
                p1.SizeData                 = 300;
                p1.MarkerFaceColor          = color;
                p1.LineWidth                = 2;
                p1.MarkerEdgeColor          = 'k';
            elseif kk==2
                p(patient) = scatter(MRD_mean,T_relapse_mean,'filled','DisplayName',text_legend_name); hold on;
                if patient<=6
                    p(patient).Marker   = 's';
                else
                    p(patient).Marker   = 'o';
                end
                p(patient).SizeData         = 300;
                p(patient).MarkerFaceColor  = color;
                p(patient).LineWidth         = 2;
                p(patient).MarkerEdgeColor   = 'k';
            else
                p1 = scatter(MRD_mean,T_relapse_mean,'>','filled'); hold on;
                p1.SizeData                 = 300;
                p1.MarkerFaceColor          = color;
                p1.LineWidth                = 2;
                p1.MarkerEdgeColor          = 'k';
            end
%           Save the results
            if (patient<=N_patients_fit)
                N_count_mean    = N_count_mean+1;
                MRD_mean_vec(N_count_mean)          = MRD_mean;
                T_relapse_mean_vec(N_count_mean)    = T_relapse_mean;
            end
        end
    end

    [MRD_values,T_relapse_values]=Hill_fit;
    p(N_patients+1) = plot(MRD_values,T_relapse_values,'-r','DisplayName','Fitted Hill function using Ding dataset'); hold on
    p(N_patients+1).LineWidth = 2;

    MRD_mean_vec(3*(i_exclude-1)+1:3*i_exclude)       = [];
    T_relapse_mean_vec(3*(i_exclude-1)+1:3*i_exclude) = [];
    [MRD_values,T_relapse_values]=Hill_fit;
    p(N_patients+2) = plot(MRD_values,T_relapse_values,'-k','DisplayName',['Fitted Hill function using Ding dataset without patient ' num2str(i_exclude)]); hold on
    p(N_patients+2).LineWidth = 2;

    set(gca,'XScale','log');
    set(gca,'fontsize',20);
    l = xlabel(['MRD at ' num2str(L_MRD_wait/30) ' months after remission (%)']);
    l.FontSize  = 25;
    l   = ylabel('Time to relapse (days)');
    l.FontSize  = 25;

    [l,objh]    = legend(p(1:N_patients+2),'Location','northeast');
    l.FontSize  = 20;
    ch = findobj(objh,'type','patch');
    set(ch,'Markersize',20);

    figure_name = ['plot_MRD_time_to_relapse_all.png'];
    saveas(gcf,figure_name);

end
%=======================================================================
%=======================================================================
%=======================================================================
function [MRD_values,T_relapse_values] = Hill_fit()
    global MRD_mean_vec T_relapse_mean_vec

    MRD_mean_vec(isnan(MRD_mean_vec))=10;
    T_relapse_mean_vec(isnan(T_relapse_mean_vec))=100;

    fun         = @(parameter_set,MRD) Hill_func(parameter_set,MRD);
    [para,resnorm]= lsqcurvefit(fun,[1600,1,5],MRD_mean_vec,T_relapse_mean_vec);
    % ,'MaxIterations',1000
    % para    = real(para);
    MRD_values  = 10.^[-7:0.1:2];
    A           = para(1);
    B           = para(2);
    C           = -7;
    n           = para(3);
    T_relapse_values= A./(B+(abs(log10(MRD_values)-C)).^n);
    fprintf('A=%f\nB=%f\nn=%f\n',para(1),para(2),para(3));
end
%=======================================================================
%=======================================================================
%=======================================================================
function T_relapse_output = Hill_func(parameter_set,MRD)
    A   = parameter_set(1);
    B   = parameter_set(2);
    C   = -7;
    n   = parameter_set(3);
    T_relapse_output    = A./(B+(abs(log10(MRD)-C)).^n);
    T_relapse_output(isnan(T_relapse_output))=0.03;
end
