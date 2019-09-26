figure(1)
clf
%--------------------------------------------Parameters for each patient
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
%-----------------------------------------------------------------------
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
%-----------------------------------------------------------------------
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
%-----------------------------------------------------------------------
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
%-----------------------------------------------------------------------
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
%=======================================================================
%   Maximum outcomes to plot
    max_outcome     = 4;
%   Number of trajectories per parameter set
    N_traj          = 1000;
%   Number of parameter sets
    NO_PARAMETERS   = 1;
%   Threshold between important and unimportant clones (only mitotic
%   cells concerned)
    N_thres         = 10^3;
%   Parameters of hematopoietic clone:
    a_c             = 0.87;
    p_c             = 0.45;
    d_c             = 2.3;
%   Equilibrium hematopoietic mitotic cell numbers:
    c_1_bar         = 2*(10^9)*weight;
%   Total number of cells in bone marrow at primary/diagnosis:
    full_BM_sum     = 4.6*(10^11);
    diag_BM_sum     = p_cellularity*full_BM_sum;
%   Constant for negative feedback of hematopoietic clone, as in
%   feedback = 1/(1+k*c2) where c2 is number of post-mitotic hematopoietic
%   cells number:
    k               = 10^-12;
%   Constants for negative feedback of overcrowding bone marrow, as in
%   feedback = C1*max(0,x-C2*c_1_d_bar) where x is total number of cells in
%   bone marrow (including mitotic hematopoietic and mitotic and 1st stage
%   post-mitotic cancer cells):
    C2              = 1;
    C1              = 10^-12;
%   Migration rate of post-mitotic cancer cells to peripheral blood:
    m               = 10;
%=====================================PLOTTING THE DETERMINISTIC RESULTS
    for i_para=1:3
%-------Input the clonal parameters
        filename        = ['parameters_' patient_ID];
        if i_para==1
            filename    = [filename '_left.txt'];
        elseif i_para==2
            filename    = [filename '_centered.txt'];
        elseif i_para==3
            filename    = [filename '_right.txt'];
        end
        fileID          = fopen(filename,'r');
        formatSpec      = '%f';
        A               = fscanf(fileID,formatSpec);
        fclose(fileID);
        A               = reshape(A,2*N_clones,[])';
        para_a_l        = A(index_parameter,1:N_clones);
        para_p_l        = A(index_parameter,N_clones+1:2*N_clones);
        para_d_l        = 0.5*ones(1,length(para_a_l));
%-------Set up parameters at diagnosis
        parameters_l= [para_a_l;  para_p_l;  para_d_l];
%  	    Population vector at diagnosis
        P           = zeros(1,2*N_clones+2);
        P(1)        = diag_BM_sum*diag_clonal(end)/100;
        P(2)        = diag_normal;
        for j=1:N_clones
%           Compute number of mitotic cells in BM
            P(2*j+1)= diag_BM_sum*diag_clonal(j)/100;
%           Compute number of mature cells in blood
            P(2*j+2)= diag_blast*diag_clonal(j)/100;
        end
%-------From diagnosis to end of chemotherapy
%       Run the ODEs for induction
        T_temp  = 0;
        ODE     = @(t,P)ODE_chemo(t,P,parameters_l,m,k,C1,C2,...
                                      k_induction,c_1_bar,a_c,p_c,d_c);
        [T_result,P_result]=ode23(ODE,[T_temp (T_temp+L_induction)],P);
        T_temp  = T_temp+L_induction;
        P_induc = P_result(end,:);
        P_temp  = P_result(end,:);
%       Save the result from the ODEs
        T_entire= T_result;
        P_entire= P_result;
%       Run the ODEs for the gap between induction and consolidation
        ODE     = @(t,P)ODE_normal(t,P,parameters_l,m,k,C1,C2,...
                                   c_1_bar,a_c,p_c,d_c);
        [T_result,P_result]=ode23(ODE,[T_temp (T_temp+30)],P_temp);
        T_temp  = T_temp+30;
        P_biopsy= P_result(end,:);
        P_temp  = P_result(end,:);
%       Save the result from the ODEs
        T_entire= [T_entire' T_result']';
        P_entire= [P_entire' P_result']';
%       Calculate clonal and normal percentages in BM
        percent_BM_Biopsy    = zeros(N_clones+1,1);
        Sum_BM  = sum(P_biopsy(2*i+1:2:end))+P_biopsy(1);
        for j=1:N_clones
            percent_BM_Biopsy(j)=100*P_biopsy(2*j+1)/Sum_BM;
        end
        percent_BM_Biopsy(N_clones+1)=100*P_biopsy(1)/Sum_BM;
%       Run the ODEs for consolidation
        ODE     = @(t,P)ODE_chemo(t,P,parameters_l,m,k,C1,C2,...
                                  k_consolidation,c_1_bar,a_c,p_c,d_c);
        [T_result,P_result]=ode23(ODE,[T_temp (T_temp+L_consolidation)],P_temp);
        T_temp  = T_temp+L_consolidation;
        P_conso = P_result(end,:);
        P_temp  = P_result(end,:);
%       Save the result from the ODEs
        T_entire= [T_entire' T_result']';
        P_entire= [P_entire' P_result']';
%       Calculate clonal and normal percentages in BM
        percent_BM_Conso     = zeros(N_clones+1,1);
        Sum_BM  = sum(P_conso(3:2:end))+P_conso(1);
        for j=1:N_clones
            percent_BM_Conso(j)=100*P_conso(2*j+1)/Sum_BM;
        end
        percent_BM_Conso(N_clones+1)=100*P_conso(1)/Sum_BM;
%-------From end of chemotherapy to relapse
%       Run the ODEs
        ODE     = @(t,P)ODE_normal(t,P,parameters_l,m,k,C1,C2,...
                                   c_1_bar,a_c,p_c,d_c);
        [T_result,P_result]=ode23(ODE,[T_temp L_relapse],P_temp);
        P_relap = P_result(end,:);
%       Save the result from the ODEs
        T_entire= [T_entire' T_result']';
        P_entire= [P_entire' P_result']';
%       Calculate clonal and normal percentages in BM
        percent_BM_Relap     = zeros(N_clones+1,1);
        Sum_BM  = sum(P_relap(3:2:end))+P_relap(1);
        for j=1:N_clones
            percent_BM_Relap(j)=100*P_relap(2*j+1)/Sum_BM;
        end
        percent_BM_Relap(N_clones+1)=100*P_relap(1)/Sum_BM;
%-------Find the bone marrow clonal percentages and cellularity over time
        BM_cellular = zeros(1,length(T_entire));
        BM_clonal   = zeros(length(T_entire),N_clones+1);
        for i=1:length(T_entire)
            BM_cellular(i)  = P_entire(i,1)+sum(P_entire(i,3:2:end));
            for j=1:N_clones
                BM_clonal(i,j) = 100*P_entire(i,2*j+1)/BM_cellular(i);
            end
            BM_clonal(i,N_clones+1)  = 100*P_entire(i,1)/BM_cellular(i);
            BM_cellular(i)  = 100*BM_cellular(i)/full_BM_sum;
        end
        max_BM_cellular = 20*ceil(max(BM_cellular)/20);
%-------Plot the clonal percentages over time
        plot_start  = (i_para-1)*20+3;
        plot_end    = i_para*20-2;
        subplot(max_outcome+1,60,[plot_start:plot_end]);
        p = area(T_entire,BM_clonal);
        set(p,'XData',T_entire);
        p(N_clones+1).FaceColor     = 'r';
        for j=1:N_clones
            if j==1     p(j).FaceColor  = 'b';
            elseif j==2 p(j).FaceColor  = 'm';
            elseif j==3 p(j).FaceColor  = 'k';
            elseif j==4 p(j).FaceColor  = 'c';
            elseif j==5 p(j).FaceColor  = [0.85 0.85 0.85];
            end
        end
        set(gca,'ytick',[])
        % set(gca,'xticklabel',[]);
        set(gca,'xaxisLocation','top');

        xlim([0 L_relapse]);    ylim([0 100]);
        a   = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',12);
        t = title(['Parameter set ' num2str(i_para)]);
        t.FontSize  = 15;
%-------Plot the diagnosis clonal percentages
        subplot(max_outcome+1,60,plot_start-1);
        b = bar([diag_clonal';zeros(1,N_clones+1)],'stacked');
        xlim([0.75 1.25]);
        b(N_clones+1).FaceColor = 'r';
        for j=1:N_clones
            if j==1     b(j).FaceColor  = 'b';
            elseif j==2 b(j).FaceColor  = 'm';
            elseif j==3 b(j).FaceColor  = 'k';
            elseif j==4 b(j).FaceColor  = 'c';
            elseif j==5 b(j).FaceColor  = [0.85 0.85 0.85];
            end
        end
        set(gca,'xtick',[]);
        a   = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',12);
%-------Plot the relapse clonal percentages
        subplot(max_outcome+1,60,plot_end+1);
        b = bar([rel_clonal';zeros(1,N_clones+1)],'stacked');

        set(gca,'XScale','log');

        xlim([0.75 1.25]);
        b(N_clones+1).FaceColor = 'r';
        for j=1:N_clones
            if j==1      b(j).FaceColor  = 'b';
            elseif j==2 b(j).FaceColor  = 'm';
            elseif j==3 b(j).FaceColor  = 'k';
            elseif j==4 b(j).FaceColor  = 'c';
            elseif j==5 b(j).FaceColor  = [0.85 0.85 0.85];
            end
        end
        ax = gca;
        ax.YAxisLocation = 'right';
        set(gca,'xtick',[]);
        a   = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',12);
        set(gca,'YTick',[]);
    end
%========================================PLOTTING THE STOCHASTIC RESULTS
    for i_para=1:3
%-------Parameters of the algorithm
        plot_1_max      = 0;
        plot_2_max      = 0;
        plot_1_min      = Inf;
        plot_2_min      = Inf;
        logic_start     = 0;
        vec_cured       = zeros(1,N_clones+1);
        vec_cured(end)  = 1;
%-------Preprocessing the outcomes
        N_count         = 0;
        % for jj=1:NO_PARAMETERS
            for ii=1:N_traj
%               Download the trajectory
                filename    = [patient_ID '/'];
                filename    = [filename '/stochastic_' patient_ID '_' num2str(i_para) '_' num2str(ii) '.txt'];
                if exist(filename,'file')~=2
                    continue;
                end
                fprintf('Preprocessing %d---%d---%d\n',i_para,ii,i_para);
                fileID      = fopen(filename,'r');
                formatSpec  = '%f';
                A           = fscanf(fileID,formatSpec);
                fclose(fileID);
                A           = reshape(A,1+2+2*N_clones,[])';
                T_entire    = A(:,1);
                if T_entire(end)~=L_relapse
                    continue;
                end
                P_entire    = A(:,2:end);
                N_count     = N_count+1;
%               Replace 0's in P_entire with 1's
                P_entire(find(P_entire<1))=1;
%               Find the important and unimportant clones at the end
                P_end       = P_entire(end,:);
                vec_char    = zeros(1,N_clones+1);
                for j=1:N_clones
                    if P_end(2*j+1)>N_thres
                        vec_char(j) = 1;
                    end
                end
                if P_end(1)>N_thres
                    vec_char(N_clones+1) = 1;
                end
%               Find which figure the trajectory belongs to
                if logic_start==0
%                   This is the first trajectory to be plotted
                    logic_start     = 1;
                    Table_char      = vec_char;
                    count_vec       = [0];
                end
                i_figure    = 0;
                for j=1:size(Table_char,1)
                    if isequal(vec_char,Table_char(j,:))
                        i_figure    = j;
                        count_vec(j)= count_vec(j)+1;
                        break;
                    end
                end
                if i_figure==0
                    Table_char(end+1,:) = vec_char;
                    count_vec(end+1)    = 1;
                    if isequal(vec_char,vec_cured)
%                       If there's a cured outcome, put it in 1st position
                        T_char(1,:)                     = Table_char(end,:);
                        c_vec(1,:)                      = count_vec(end);
                        T_char(2:length(count_vec),:)   = Table_char(1:end-1,:);
                        c_vec(2:length(count_vec))      = count_vec(1:end-1);
                        Table_char                      = T_char;
                        count_vec                       = c_vec;
                    end
                end
            end
        % end
        disp(N_count);
%-------Plotting the outcomes
        N_outcome   = min(length(count_vec),max_outcome);
        % for jj=1:NO_PARAMETERS
            for ii=1:N_traj
%               Download the trajectory
                filename    = [patient_ID '/'];
                filename    = [filename '/stochastic_' patient_ID '_' num2str(i_para) '_' num2str(ii) '.txt'];
                if exist(filename,'file')~=2
                    continue;
                end
                fprintf('Plotting %d---%d---%d\n',i_para,ii,i_para);
                fileID      = fopen(filename,'r');
                formatSpec  = '%f';
                A           = fscanf(fileID,formatSpec);
                fclose(fileID);
                A           = reshape(A,1+2+2*N_clones,[])';
                T_entire    = A(:,1);
                if T_entire(end)~=L_relapse
                    continue;
                end
                P_entire    = A(:,2:end);
%               Replace 0's in P_entire with 1's
                P_entire(find(P_entire<1))=1;
%               Find the important and unimportant clones at the end
                P_end       = P_entire(end,:);
                vec_char    = zeros(1,N_clones+1);
                for j=1:N_clones
                    if P_end(2*j+1)>N_thres
                        vec_char(j) = 1;
                    end
                end
                if P_end(1)>N_thres
                    vec_char(N_clones+1) = 1;
                end
%               Find which figure the trajectory belongs to
                i_figure    = 0;
                for j=1:size(Table_char,1)
                    if isequal(vec_char,Table_char(j,:))
                        i_figure    = j;
                        break;
                    end
                end
                if i_figure>N_outcome
                    continue;
                end
%               Some parameters for the plots
                plot_1_max  = max(plot_1_max,1.2*max(max(P_entire(:,1)),max(max(P_entire(:,3:2:end)))));
                plot_1_min  = min(plot_1_min,1.2*min(min(P_entire(:,1)),min(min(P_entire(:,3:2:end)))));
%               Plot the log-plot in BM
                plot_start  = i_figure*60+(i_para-1)*20+3;
                plot_end    = i_figure*60+i_para*20-2;
                subplot(max_outcome+1,60,[plot_start:plot_end]);
                p = semilogy(T_entire,P_entire(:,1)); hold on
                p.LineWidth     = 2;
                p.Color         = 'r';
                for j=1:N_clones
                    p = semilogy(T_entire,P_entire(:,2*j+1)); hold on
%                     T_entire

                    set(gca,'XScale','log')

                    p.LineWidth     = 2;
                    if j==1         p.Color     = 'b';
                        elseif j==2 p.Color     = 'm';
                        elseif j==3 p.Color     = 'k';
                        elseif j==4 p.Color     = 'c';
                        elseif j==5 p.Color     = [0.85 0.85 0.85];
                    end
                end
            end
        % end
%-------Illustrating the plots
        for j=1:N_outcome
            plot_start  = j*60+(i_para-1)*20+3;
            plot_end    = j*60+i_para*20-2;
            subplot(max_outcome+1,60,[plot_start:plot_end]);

            % subplot(4,40,[(40*j-16):(40*j)]);
            p = rectangle('Position',[1 1 (L_induction-1) 10^ceil(log10(plot_1_max))],'FaceColor','g');hold on
                uistack(p,'bottom');
            p = rectangle('Position',[(L_induction+30) 1 L_consolidation 10^ceil(log10(plot_1_max))],'FaceColor','g');hold on
                uistack(p,'bottom');
            xlim([0 L_relapse]);    ylim([max(1,10^floor(log10(plot_1_min))) 10^ceil(log10(plot_1_max))]);
            vecYTick = [max(1,floor(log10(plot_1_min))):3:ceil(log10(plot_1_max))];
            for i=1:length(vecYTick)
                vecYTick(i) = 10^vecYTick(i);
            end
            ax = gca;
            ax.YTick = vecYTick;
            perc_sto    = 100*count_vec(j)/N_count;

            fprintf('\nCase %d with probability %f:\n',j,perc_sto);
            disp(Table_char(j,:));
            fprintf('\n');

            vec_alive   = find(Table_char(j,:)==1);
            N_alive     = length(vec_alive);
            vec_dead    = find(Table_char(j,:)==0);
            N_dead      = N_clones+1-N_alive;

            a   = get(gca,'XTickLabel');
            set(gca,'XTickLabel',a,'fontsize',12);
            Title   = [num2str(round(perc_sto,1)) '%'];
            t = title(Title);
            t.FontSize  = 15;
            if j==N_outcome
                xlabel('days');
            else
                set(gca,'xticklabel',[]);
            end
            ylabel('cells');
        end


    end

    figure_name = ['plot_stochastic_' patient_ID '.png'];
    saveas(gcf,figure_name);

%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
function dP=ODE_normal(t,P,parameters_l,m,k,C1,C2,c_1_bar,a_c,p_c,d_c)
%   ABSTRACT----This function models the ODE during the normal time (not
%               under chemotherapy).
%   IMPORTANT---* The part where mitotic populations falling under 1 cell
%               are killed might need to be changed to be more efficient.
%   INPUT-------t               [REAL] input time to compute the ODEs.
%               P               [VECTOR OF LENGTH (3*n+2) OF POSITIVE REAL]
%                               population vector at the moment, where n is
%                               the number of leukemic clones.
%               parameters_l    [MATRIX OF SIZE (3)x(n) of POSITIVE REAL]
%                               parameters for the leukemic clones. First
%                               row is renewal rates (<=1), second row is
%                               proliferation rates, last row is death
%                               rates.
%               m               [POSITIVE REAL] migration rate of mature
%                               leukemic cells from bone marrow to blood.
%               k               [POSITIVE REAL] feedback constant for
%                               overpopulation in blood.
%               C1,C2           [POSITIVE REAL] feedback constants for
%                               overpopulation in bone marrow.
%               c_1_bar         [POSITIVE REAL] number of normal mitotic
%                               cells at equilibrium.
%               a_c,p_c,d_c     [POSITIVE REAL] parameters for the normal
%                               clone. Note that a_c<1.
%   OUTPUT------dP              [VECTOR OF LENGTH (3*n+2)] 1st derivatives
%                               of cells numbers according to Model 3.
%   LAST CHECK--June 13, 2016.
%-----------------------------------------------Parameters of the algorithm
%   Number of clones
    N_clones    = size(parameters_l,2);
%--------------------------------------Parameters in the mathematical model
    Sum         = sum(P(4:2:2*N_clones+2))+P(2);
    s           = 1/(1+k*Sum);
    x           = sum(P(3:2:2*N_clones+2))+P(1);
    d           = C1*max(0,x-C2*c_1_bar);
%-------------------------------------------------------Deterministic model
%   Hematopoietic clone: (if mitotic cells fall under 1, they're dead)
    dP(1) = (2*a_c*s-1)*p_c*P(1)-d*P(1);
    dP(2) = 2*(1-a_c*s)*p_c*P(1)-d_c*P(2);
    if P(1)<0.95
        dP(1)=0;
    end
%   Leukemic clones: (if mitotic cells fall under 1, they're dead)
    for i=1:N_clones
        a_l = parameters_l(1,i);
        p_l = parameters_l(2,i);
        d_l = parameters_l(3,i);
        dP(2*i+1) = (2*a_l*s-1)*p_l*P(2*i+1)-d*P(2*i+1);
        dP(2*i+2) = 2*(1-a_l*s)*p_l*P(2*i+1)-d_l*P(2*i+2);
        if P(2*i+1)<0.8
            dP(2*i+1)=0;
        end
    end
%------------------------------------------------------Transpose the vector
    dP=dP';
end
%==========================================================================
function dP=ODE_chemo(t,P,parameters_l,m,k,C1,C2,...
                          k_chemo,c_1_bar,a_c,p_c,d_c)
%   ABSTRACT----This function models the ODE during the treatment time
%               (under chemotherapy).
%   IMPORTANT---* The part where mitotic populations falling under 1 cell
%               are killed might need to be changed to be more efficient.
%   INPUT-------t               [REAL] input time to compute the ODEs.
%               P               [VECTOR OF LENGTH (3*n+2) OF POSITIVE REAL]
%                               population vector at the moment, where n is
%                               the number of leukemic clones.
%               parameters_l    [MATRIX OF SIZE (3)x(n) of POSITIVE REAL]
%                               parameters for the leukemic clones. First
%                               row is renewal rates (<=1), second row is
%                               proliferation rates, last row is death
%                               rates.
%               m               [POSITIVE REAL] migration rate of mature
%                               leukemic cells from bone marrow to blood.
%               k               [POSITIVE REAL] feedback constant for
%                               overpopulation in blood.
%               C1,C2           [POSITIVE REAL] feedback constants for
%                               overpopulation in bone marrow.
%               k_chemo         [VECTOR OF LENGTH 2 OF POSITIVE REAL]
%                               chemotherapy constants. The first number is
%                               for normal clone, the second one is for
%                               leukemic clones.
%               c_1_bar         [POSITIVE REAL] number of normal mitotic
%                               cells at equilibrium.
%               a_c,p_c,d_c     [POSITIVE REAL] parameters for the normal
%                               clone. Note that a_c<1.
%   OUTPUT------dP              [VECTOR OF LENGTH (3*n+2)] 1st derivatives
%                               of cells numbers according to Model 3.
%   LAST CHECK--June 13, 2016.
%-----------------------------------------------Parameters of the algorithm
%   Number of clones
    N_clones    = size(parameters_l,2);
%---------------------------------------------Get result from normal period
    dP = ODE_normal(t,P,parameters_l,m,k,C1,C2,c_1_bar,a_c,p_c,d_c);
%----------------------------------Additional death rates for mitotic cells
    dP(1)       = dP(1)-k_chemo(1)*p_c*P(1);
    for i=1:N_clones
        p_l = parameters_l(2,i);
        dP(2*i+1) = dP(2*i+1)-k_chemo(2)*p_l*P(2*i+1);
    end
end
