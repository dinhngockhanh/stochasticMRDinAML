function PROJECT_growth_rate
%-------------------------------------------Parameters for the algorithm
%   Number of patients
    N_patients      = 6;
%   Weights of all patients
    weight          = [70.3   87.5   103.6  80.7   90.2   79.3];
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
    p_c                 = zeros(max(N_cancer),6);
    a_c                 = zeros(max(N_cancer),6);
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
        % para_d_l        = 0.2*ones(1,length(para_a_l));
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
%   Constants for two negative feedback systems:
    C1                  = 10^-12;
    C2                  = 1;
    k                   = 10^-12;
%   Plot background
    figure(1)
    clf
    plot([0 6],[0.3 0.3],'--r','LineWidth',2);hold on
    plot([0 6],[2.0 2.0],'--r','LineWidth',2);hold on
%------------------------------------------------------------Main algorithm
    for i=1:6
        % fprintf('Patient %d:\n',i);
        c_1_bar = 2*(10^9)*weight(i);
        s_bar   = 1/(1+k*normal_in_blood(i));
        d_bar   = C1*(normal_in_BM(i)-C2*c_1_bar);
        for j=1:N_cancer(i)
            lambda = 30*(log10(exp(1)))*((2*a_c(j,i)*s_bar-1)*p_c(j,i)-d_bar)
            scatter([i],[lambda],[100],'MarkerEdgeColor',[0 .5 .5],...
                                 'MarkerFaceColor',[0 .7 .7],...
                                 'LineWidth',3);hold on
                             % ylim([0 2.5]);
                             xlim([0 6]);
        end
    end
    l = xlabel('Patient');
    l.FontSize  = 15;
    l = ylabel('Growth rate (lambda)');
    l.FontSize  = 15;

    figure_name = ['plot_growth_rate_all.png'];
    saveas(gcf,figure_name);
end
