function matlab_plotting_MRD_analysis
%-----------------------------------------------Parameters for each patient
% patient_ID          = '400220';
% index_parameter     = 1;
% N_clones            = 2;
% weight              = 70.3;
% p_cellularity       = 0.9;
% diag_normal         = 2.94*(10^9);
% diag_blast          = 5.68*(10^10);
% diag_clonal         = [68   3   29]';
% rel_clonal          = [38   22  40]';
% L_induction         = 7;
% k_induction         = [1*0.8 5*0.8];
% L_consolidation     = 11;
% k_consolidation     = [1*0.8 5*0.8];
% L_relapse           = 252;
% BMblast_relapse     = 60;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% patient_ID          = '426980';
% index_parameter     = 1;
% N_clones            = 5;
% weight              = 87.5;
% p_cellularity       = 0.5;
% diag_normal         = 1.65*(10^9);
% diag_blast          = 1.32*(10^10);
% diag_clonal         = [3    3   33  22  3   36]';
% rel_clonal          = [0    3   0   0   9   88]';
% L_induction         = 7;
% k_induction         = [1*0.9 5*0.9];
% L_consolidation     = 8;
% k_consolidation     = [1*0.9 5*0.9];
% L_relapse           = 805;
% BMblast_relapse     = 12;
%--------------------------------------------------------------------------
% patient_ID          = '452198';
% index_parameter     = 1;
% N_clones            = 4;
% weight              = 103.6;
% p_cellularity       = 0.9;
% diag_normal         = 3.63*(10^9);
% diag_blast          = 2.90*(10^10);
% diag_clonal         = [19   52  23  3   3]';
% rel_clonal          = [0    0   0   20  80]';
% L_induction         = 7;
% k_induction         = [1*1.0 5*1.0];
% L_consolidation     = 11;
% k_consolidation     = [1*1.0 5*1.0];
% L_relapse           = 505;
% BMblast_relapse     = 20;
%--------------------------------------------------------------------------
% patient_ID          = '573988';
% index_parameter     = 1;
% N_clones            = 2;
% weight              = 80.7;
% p_cellularity       = 0.9;
% diag_normal         = 1.22*(10^10);
% diag_blast          = 7.60*(10^9);
% diag_clonal         = [72   3   25]';
% rel_clonal          = [15   39  46]';
% L_induction         = 7;
% k_induction         = [1*1.0 5*1.0];
% L_consolidation     = 5;
% k_consolidation     = [1*1.0 5*1.0];
% L_relapse           = 365;
% BMblast_relapse     = 54;
%--------------------------------------------------------------------------
% patient_ID          = '758168';
% index_parameter     = 1;
% N_clones            = 3;
% weight              = 90.2;
% p_cellularity       = 0.9;
% diag_normal         = 3.26*(10^9);
% diag_blast          = 0.0;
% diag_clonal         = [7    83  3   7]';
% rel_clonal          = [39   18  35  8]';
% L_induction         = 21;
% k_induction         = [1*0.15 5*0.15];
% L_consolidation     = 90;
% k_consolidation     = [1*0.15 5*0.15];
% L_relapse           = 961;
% BMblast_relapse     = 92;
%--------------------------------------------------------------------------
patient_ID          = '804168';
index_parameter     = 1;
N_clones            = 2;
weight              = 79.3;
p_cellularity       = 0.9;
diag_normal         = 1.54*(10^11);
diag_blast          = 2.29*(10^11);
diag_clonal         = [83   3   14]';
rel_clonal          = [36   45  19]';
L_induction         = 7;
k_induction         = [1*0.8 5*0.8];
L_consolidation     = 11;
k_consolidation     = [1*0.8 5*0.8];
L_relapse           = 235;
BMblast_relapse     = 81;
%=======================================================================
%   Time of remission
    L_remission         = L_induction+30+L_consolidation;
%   Number of trajectories per parameter set
    N_traj              = 1000;
%   Number of parameter sets
    NO_PARAMETERS       = 3;
%   Level of "detectible" MRD
    level_detect_MRD    = 0.01;
%=======================================================================
%   Prepare the arrays for plotting
    vec_time            = [0:1:L_relapse];

    vec_PBblast_mean    = zeros(1,length(vec_time));
    vec_PBblast_min     = 100*ones(1,length(vec_time));
    vec_PBblast_max     = zeros(1,length(vec_time));

    vec_BMblast_mean    = zeros(1,length(vec_time));
    vec_BMblast_min     = 100*ones(1,length(vec_time));
    vec_BMblast_max     = zeros(1,length(vec_time));

    vec_time_5perc      = [0:10:3*L_relapse];
    T_relapse           = zeros(1,length(vec_time_5perc));

    MRD_relapse_scale   = [0:1:100];
    MRD_relapse         = zeros(1,length(MRD_relapse_scale));

    vec_detect_MRD      = zeros(1,N_traj*NO_PARAMETERS);
%--------------------
    N_count         = 0;
    N_count_2       = 0;
    index_detectible_MRD    = 0;
    for jj=1:NO_PARAMETERS
        for ii=1:N_traj
%           Download the trajectory
            filename    = [patient_ID '/stochastic_' patient_ID '_' num2str(jj) '_' num2str(ii) '.txt'];
            % filename    = [patient_ID '/stochastic_' patient_ID '_' num2str(jj) '_' num2str(ii) '.txt'];
            if exist(filename,'file')~=2
                continue;
            end
            fprintf('Processing result %d of parameter %d\n',ii,jj);
            fileID      = fopen(filename,'r');
            formatSpec  = '%f';
            A           = fscanf(fileID,formatSpec);
            fclose(fileID);
            A           = reshape(A,1+2+2*N_clones,[])';
%           Output the time and trajectory arrays
            T_entire    = A(:,1);
            if T_entire(end)<L_relapse
                continue;
            end
            P_entire    = A(:,2:end);
%-----------Find evolution of MRD in blood and BM
            N_count = N_count+1;
            for kk=1:length(vec_time)
                index   = find(T_entire>=vec_time(kk),1);
%               Find %BM blast
                BM_blast    = 0;
                for i_clone=1:N_clones
                    BM_blast    = BM_blast+P_entire(index,2*i_clone+1);
                end
                Sum     = BM_blast+P_entire(index,1);
                BMblast = 100*BM_blast/Sum;
                vec_BMblast_mean(kk)    = vec_BMblast_mean(kk)+BMblast;
                vec_BMblast_min(kk)     = max(10^-10,min(vec_BMblast_min(kk),BMblast));
                vec_BMblast_max(kk)     = max(vec_BMblast_max(kk),BMblast);
%               Find %PB blast
                PB_blast    = 0;
                for i_clone=1:N_clones
                    PB_blast    = PB_blast+P_entire(index,2*i_clone+2);
                end
                Sum     = PB_blast+P_entire(index,2);
                PBblast = 100*PB_blast/Sum;
                vec_PBblast_mean(kk)    = vec_PBblast_mean(kk)+PBblast;
                vec_PBblast_min(kk)     = max(10^-10,min(vec_PBblast_min(kk),PBblast));
                % vec_PBblast_min(kk)     = min(vec_PBblast_min(kk),PBblast);
                vec_PBblast_max(kk)     = max(vec_PBblast_max(kk),PBblast);
            end
%-----------Find "time to relapse", defined as when %PB blast>=5% and %BM
%-----------blast at at time; and
%-----------"time to detectible MRD", defined as when %PB blast>=0.01%
            found_T_relapse         = 0;
            found_detectible_MRD    = 0;
            index                   = 0;
            while (found_T_relapse==0)
                index   = index+1;
%               Find %BM blast
                BM_blast    = 0;
                for i_clone=1:N_clones
                    BM_blast    = BM_blast+P_entire(index,2*i_clone+1);
                end
                Sum     = BM_blast+P_entire(index,1);
                BMblast = 100*BM_blast/Sum;
%               Find %PB blast
                PB_blast    = 0;
                for i_clone=1:N_clones
                    PB_blast    = PB_blast+P_entire(index,2*i_clone+2);
                end
                Sum     = PB_blast+P_entire(index,2);
                PBblast = 100*PB_blast/Sum;
%               Update if "time to detectible MRD" is found...
                if ((T_entire(index)>L_remission)&&(PBblast>level_detect_MRD)&&(found_detectible_MRD==0))
                    found_detectible_MRD = 1;
                    index_detectible_MRD = index_detectible_MRD+1;
                    vec_detect_MRD(index_detectible_MRD) = T_entire(index)-L_remission;
                end
%               Update if "time to relapse" is found...
                if ((T_entire(index)>L_remission)&&(PBblast>5))
                    found_T_relapse = 1;
                    N_count_2       = N_count_2+1;
                    index2  = find(vec_time_5perc<T_entire(index),1,'last');
                    T_relapse(index2)   = T_relapse(index2)+1;
                    index2  = find(MRD_relapse_scale<BMblast,1,'last');
                    MRD_relapse(index2) = MRD_relapse(index2)+1;
%               In case the simulation never reaches relapse...
                elseif (index==length(T_entire))
                    found_T_relapse = 1;
                end
            end
        end
    end
    fprintf('\n');
    disp(patient_ID);
    fprintf('\n1st quantile, median, and 3rd quantile\n');
    fprintf('\nof time from remission to when %%MRD=%f:   %f %f %f\n',level_detect_MRD,quantile(vec_detect_MRD,[0.25 0.50 0.75]));
%     quantile(vec_detect_MRD,[0.25 0.50 0.75])
    disp(sort(vec_detect_MRD));

    vec_BMblast_mean    = vec_BMblast_mean/N_count;
    vec_PBblast_mean    = vec_PBblast_mean/N_count;

    figure(1);
    clf
    subplot(2,4,[1:3,5:7]);

    p1 = semilogy(vec_time,vec_BMblast_mean);hold on;
        p1.Color     = 'red';
        p1.LineStyle = '-';
        p1.LineWidth = 3;
    p1m = semilogy(vec_time,vec_BMblast_min);hold on;
        p1m.Color     = 'red';
        p1m.LineStyle = '--';
        p1m.LineWidth = 1.5;
    p = semilogy(vec_time,vec_BMblast_max);hold on;
        p.Color     = 'red';
        p.LineStyle = '--';
        p.LineWidth = 1.5;
    p2 = semilogy(vec_time,vec_PBblast_mean);hold on;
        p2.Color     = 'blue';
        p2.LineStyle = '-';
        p2.LineWidth = 3;
    p2m = semilogy(vec_time,vec_PBblast_min);hold on;
        p2m.Color     = 'blue';
        p2m.LineStyle = '--';
        p2m.LineWidth = 1.5;
    p = semilogy(vec_time,vec_PBblast_max);hold on;
        p.Color     = 'blue';
        p.LineStyle = '--';
        p.LineWidth = 1.5;
    p1r = plot(L_relapse,BMblast_relapse,'ro');
        p1r.Marker          = 'o';
        p1r.MarkerEdgeColor = 'r';
        p1r.MarkerFaceColor = 'r';
        p1r.MarkerSize      = 10;
        % p1r.LineWidth       = 3;

    xlim([0 L_relapse]);
    ylim([10^-8 inf]);

    lgd = legend([p1 p1m p1r p2 p2m],{'Mean %BM blast','Min/max %BM blast',...
                                      'Clinical %BM blast',...
                                      'Mean %PB blast','Min/max %PB blast'},'Location','southeast');
    lgd.FontSize = 15;
    set(gca,'fontsize',15);
    l   = xlabel('Days');
    l.FontSize  = 20;
    l   = ylabel('Percentage');
    l.FontSize  = 20;
    t   = title('A');
    t.FontSize  = 20;

    % subplot(2,4,[4]);
    % T_relapse           = T_relapse/N_count_2;
    % index1              = find(T_relapse>0,1);
    % index2              = find(T_relapse>0,1,'last');
    % bar(vec_time_5perc(index1:index2),T_relapse(index1:index2));hold on
    % l = line([L_relapse L_relapse], get(gca, 'ylim'));hold on
    % l.LineWidth = 3;
    % l.Color     = 'r';
    % l.LineStyle = ':';

    ax  = gca;
    ax.FontSize = 15;
%   400220...
    % ax.XTick    = [0:50:260];
    % xlim([0 260]);
%   426980...
%   452198...
    % ax.XTick    = [0:250:750];
    % xlim([0 750]);
    % ylim([0 0.2]);
%   573988...
    % ax.XTick    = [0:100:400];
    % xlim([0 400]);
    % ylim([0 0.5]);
%   758168...
    % ax.XTick    = [250:250:1000];
    % xlim([0 1000]);
    % ylim([0 0.6]);
    % subplot(2,4,[1 2 3 5 6 7]); ylim([10^-1 inf]);
    % subplot(2,4,[4])
%   804168...
    % ax.XTick    = [0:50:250];
    % xlim([0 250]);
    % ylim([0 0.5]);


    % l   = xlabel(['Days']);
    % l.FontSize = 20;
    % l   = ylabel('Frequency');
    % l.FontSize = 20;
    % t   = title('B');
    % t.FontSize  = 20;

    subplot(2,4,[8]);
    MRD_relapse         = MRD_relapse/sum(MRD_relapse);
    index1              = find(MRD_relapse>0,1);
    index2              = find(MRD_relapse>0,1,'last');
    bar(MRD_relapse_scale(index1:index2),MRD_relapse(index1:index2));
    set(gca,'fontsize',15);
    l   = xlabel(['% BM blast']);
    l.FontSize = 20;
    l   = ylabel('Frequency');
    l.FontSize = 20;
    t   = title('B');
    t.FontSize  = 20;

    fprintf('\nNumber of trajectories with final MRD>5: %d\n',N_count_2);

    figure_name = ['plot_MRD_analysis_' patient_ID '.png'];
    saveas(gcf,figure_name);

end
