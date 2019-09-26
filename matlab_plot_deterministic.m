function matlab_plotting_deterministic
%-----------------------------------------------Parameters for each patient
patient_ID          = '400220';
index_parameter     = 1;
N_clones            = 2;
weight              = 70.3;
p_cellularity       = 0.9;
diag_normal         = 2.94*(10^9);
diag_blast          = 5.68*(10^10);
diag_clonal         = [68   3   29]';
rel_clonal          = [38   22  40]';
L_induction         = 7;
k_induction         = [1*0.8 5*0.8];
L_consolidation     = 11;
k_consolidation     = [1*0.8 5*0.8];
L_relapse           = 252;
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
% patient_ID          = '804168';
% index_parameter     = 1;
% N_clones            = 2;
% weight              = 79.3;
% p_cellularity       = 0.9;
% diag_normal         = 1.54*(10^11);
% diag_blast          = 2.29*(10^11);
% diag_clonal         = [83   3   14]';
% rel_clonal          = [36   45  19]';
% L_induction         = 7;
% k_induction         = [1*0.8 5*0.8];
% L_consolidation     = 11;
% k_consolidation     = [1*0.8 5*0.8];
% L_relapse           = 235;
%==========================================================================
%-------------------------------------------------Options for the algorithm
%   File name for output:
    filename        = [patient_ID '.txt'];
%   Plot or not? (1 if does, 0 if not)
    plot_or_not     = 1;
%   Parameters of hematopoietic clone:
    a_c             = 0.87;
    p_c             = 0.45;
    d_c             = 2.3;
%   The specific leukemic clones' parameters:
    filename2       = ['parameters_' patient_ID '_right.txt'];
    fileID          = fopen(filename2,'r');
    formatSpec      = '%f';
    A               = fscanf(fileID,formatSpec);
    fclose(fileID);
    A               = reshape(A,2*N_clones,[])';
    para_a_l        = A(1:N_clones);
    para_p_l        = A(N_clones+1:2*N_clones);
    para_d_l        = 0.5*ones(1,length(para_a_l));
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
%--------------------------------------------Set up parameters at diagnosis
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
%-------------------------------------From diagnosis to end of chemotherapy
%   Run the ODEs for induction
    T_temp  = 0;
    ODE     = @(t,P)ODE_chemo(t,P,parameters_l,m,k,C1,C2,...
                                  k_induction,c_1_bar,a_c,p_c,d_c);
    [T_result,P_result]=ode23(ODE,[T_temp (T_temp+L_induction)],P);
    T_temp  = T_temp+L_induction;
    P_induc = P_result(end,:);
    P_temp  = P_result(end,:);
%   Save the result from the ODEs
    T_entire= T_result;
    P_entire= P_result;
%   Run the ODEs for the gap between induction and consolidation
    ODE     = @(t,P)ODE_normal(t,P,parameters_l,m,k,C1,C2,...
                               c_1_bar,a_c,p_c,d_c);
    [T_result,P_result]=ode23(ODE,[T_temp (T_temp+30)],P_temp);
    T_temp  = T_temp+30;
    P_biopsy= P_result(end,:);
    P_temp  = P_result(end,:);
%   Save the result from the ODEs
    T_entire= [T_entire' T_result']';
    P_entire= [P_entire' P_result']';
%   Calculate clonal and normal percentages in BM
    percent_BM_biopsy    = zeros(N_clones+1,1);
    Sum_BM  = sum(P_biopsy(3:2:end))+P_biopsy(1);
    for j=1:N_clones
        percent_BM_biopsy(j)=100*P_biopsy(2*j+1)/Sum_BM;
    end
    percent_BM_biopsy(N_clones+1)=100*P_biopsy(1)/Sum_BM;
%   Run the ODEs for consolidation
    ODE     = @(t,P)ODE_chemo(t,P,parameters_l,m,k,C1,C2,...
                              k_consolidation,c_1_bar,a_c,p_c,d_c);
    [T_result,P_result]=ode23(ODE,[T_temp (T_temp+L_consolidation)],P_temp);
    T_temp  = T_temp+L_consolidation;
    P_conso = P_result(end,:);
    P_temp  = P_result(end,:);
%   Save the result from the ODEs
    T_entire= [T_entire' T_result']';
    P_entire= [P_entire' P_result']';
%   Calculate clonal and normal percentages in BM
    percent_BM_Conso     = zeros(N_clones+1,1);
    Sum_BM  = sum(P_conso(3:2:end))+P_conso(1);
    for j=1:N_clones
        percent_BM_Conso(j)=100*P_conso(2*j+1)/Sum_BM;
    end
    percent_BM_Conso(N_clones+1)=100*P_conso(1)/Sum_BM;
%---------------------------------------From end of chemotherapy to relapse
%   Run the ODEs
    ODE     = @(t,P)ODE_normal(t,P,parameters_l,m,k,C1,C2,...
                               c_1_bar,a_c,p_c,d_c);
    [T_result,P_result]=ode23(ODE,[T_temp L_relapse],P_temp);
    P_relap = P_result(end,:);
%   Save the result from the ODEs
    T_entire= [T_entire' T_result']';
    P_entire= [P_entire' P_result']';
%   Calculate clonal and normal percentages in BM
    percent_BM_Relap     = zeros(N_clones+1,1);
    Sum_BM  = sum(P_relap(3:2:end))+P_relap(1);
    for j=1:N_clones
        percent_BM_Relap(j)=100*P_relap(2*j+1)/Sum_BM;
    end
    percent_BM_Relap(N_clones+1)=100*P_relap(1)/Sum_BM;
%-------------------------------------------------Plot the clonal evolution
    if plot_or_not==1
%   Some parameters for the plots
    plot_1_max  = 1.2*max(max(P_entire(:,1)),max(max(P_entire(:,3:2:end))));
    plot_2_max  = 1.2*max(max(P_entire(:,2)),max(max(P_entire(:,4:2:end))));
    plot_1_min  = 1.2*min(min(P_entire(:,1)),min(min(P_entire(:,3:2:end))));
    plot_2_min  = 1.2*min(min(P_entire(:,2)),min(min(P_entire(:,4:2:end))));
%   Plot the clonal evolution in semilogy
    clf
    l = cell(1,N_clones+1);
    subplot(5,20,[2:19]);
        h(N_clones+1)               = semilogy(T_entire,P_entire(:,1)); hold on
        h(N_clones+1).LineWidth     = 2;
        h(N_clones+1).Color         = 'r';
        l{N_clones+1}               = sprintf('normal clone');
    subplot(5,20,[22:39]);
        p = semilogy(T_entire,P_entire(:,2)); hold on
        p.LineWidth     = 2;
        p.Color         = 'r';
    for j=1:N_clones
        subplot(5,20,[2:19]);
            h(j)    = semilogy(T_entire,P_entire(:,2*j+1)); hold on
            h(j).LineWidth     = 2;
            if j==1     h(j).Color     = 'b';
            elseif j==2 h(j).Color     = 'm';
            elseif j==3 h(j).Color     = 'k';
            elseif j==4 h(j).Color     = 'c';
            elseif j==5 h(j).Color     = [0.85 0.85 0.85];
            end
            l{j}	= sprintf(['clone ' num2str(j)]);
        subplot(5,20,[22:39]);
            p = semilogy(T_entire,P_entire(:,2*j+2)); hold on
            p.LineWidth     = 2;
            if j==1     p.Color     = 'b';
            elseif j==2 p.Color     = 'm';
            elseif j==3 p.Color     = 'k';
            elseif j==4 p.Color     = 'c';
            elseif j==5 p.Color     = [0.85 0.85 0.85];
            end
    end
    subplot(5,20,[2:19]);
        p = rectangle('Position',[0 1 L_induction 10^ceil(log10(plot_1_max))],'FaceColor','g');hold on
        uistack(p,'bottom');
        p = rectangle('Position',[(L_induction+30) 1 L_consolidation 10^ceil(log10(plot_1_max))],'FaceColor','g');hold on
        uistack(p,'bottom');
        xlim([0 L_relapse]);    ylim([max(1,10^floor(log10(plot_1_min))) 10^ceil(log10(plot_1_max))]);
        vecYTick = [floor(log10(plot_1_min)):3:ceil(log10(plot_1_max))];
        for i=1:length(vecYTick)
            vecYTick(i) = 10^vecYTick(i);
        end
        ax = gca;
        ax.YTick = vecYTick;
        a   = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',12);
        t = title('A');
        t.FontSize  = 18;
        legend(h,l,'Location','southeast');
    subplot(5,20,[22:39]);
        p = rectangle('Position',[0 1 L_induction 10^ceil(log10(plot_2_max))],'FaceColor','g');hold on
        uistack(p,'bottom');
        p = rectangle('Position',[(L_induction+30) 1 L_consolidation 10^ceil(log10(plot_2_max))],'FaceColor','g');hold on
        uistack(p,'bottom');
        xlim([0 L_relapse]);    ylim([max(1,10^floor(log10(plot_2_min))) 10^ceil(log10(plot_2_max))]);
        vecYTick = [floor(log10(plot_2_min+1)):3:ceil(log10(plot_2_max))];
        for i=1:length(vecYTick)
            vecYTick(i) = 10^vecYTick(i);
        end
        ax = gca;
        ax.YTick = vecYTick;
        a   = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',12);
        t = title('B');
        t.FontSize  = 18;
%   Find the bone marrow clonal percentages and cellularity over time
    for i=1:length(T_entire)
        BM_total        = P_entire(i,1)+sum(P_entire(i,3:2:end));
        PB_total        = P_entire(i,2)+sum(P_entire(i,4:2:end));
        for j=1:N_clones
            BM_clonal(i,j)  = 100*P_entire(i,2*j+1)/BM_total;
            PB_clonal(i,j)  = 100*P_entire(i,2*j+2)/PB_total;
        end
        BM_clonal(i,N_clones+1) = 100*P_entire(i,1)/BM_total;
        PB_clonal(i,N_clones+1) = 100*P_entire(i,2)/PB_total;

        BM_cellular(i)  = 100*BM_total/full_BM_sum;
    end
    max_BM_cellular = 20*ceil(max(BM_cellular)/20);
%   Plot the bone marrow cellularity over time
    subplot(5,20,[42:59]);
                        p = plot(T_entire,BM_cellular);
                        p.LineWidth     = 2;
                        p.Color         = [0.5 0.5 0];
                        p = rectangle('Position',[0 0 L_induction max_BM_cellular],'FaceColor','g');hold on
                        uistack(p,'bottom');
                        p = rectangle('Position',[(L_induction+30) 0 L_consolidation max_BM_cellular],'FaceColor','g');hold on
                        uistack(p,'bottom');
                        xlim([0 L_relapse]);   ylim([0 max_BM_cellular]);
                        legend('bone marrow cellularity (%)','Location','southeast');
                        vecYTick = [0:20:max_BM_cellular];

                        ax = gca;
                        ax.YTick = vecYTick;
                        a   = get(gca,'XTickLabel');
                        set(gca,'XTickLabel',a,'fontsize',12);
                        t = title('C');
                        t.FontSize  = 18;
%   Plot the clonal percentages in BM over time
    subplot(5,20,[62:79]);     p = area(BM_clonal);
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
                        xlim([0 L_relapse]);    ylim([0 100]);
                        a   = get(gca,'XTickLabel');
                        set(gca,'XTickLabel',a,'fontsize',12);
                        t = title('D');
                        t.FontSize  = 18;
%   Plot the diagnosis and relapse clonal percentages
    subplot(5,20,61); b = bar([diag_clonal';zeros(1,N_clones+1)],'stacked');
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
                       set(gca,'xtick',[]);
                       a   = get(gca,'XTickLabel');
                        set(gca,'XTickLabel',a,'fontsize',12);
    subplot(5,20,80); b = bar([rel_clonal';zeros(1,N_clones+1)],'stacked');
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
%   Plot the clonal percentages in blood over time
    subplot(5,20,[82:99]);     p = area(PB_clonal);
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
                        xlim([0 L_relapse]);    ylim([0 100]);
                        a   = get(gca,'XTickLabel');
                        set(gca,'XTickLabel',a,'fontsize',12);
                        t = title('E');
                        t.FontSize  = 18;
%   End plotting session
    figure_name = ['plot_deterministic_' patient_ID '.png'];
    saveas(gcf,figure_name);
    end
%--------------------------------------------------------Output the outcome
%   Output results on screen
    fprintf('Clonal percentages right before consolidation (biopsy):\n');
    disp(percent_BM_biopsy);
    fprintf('BM populations right after consolidation:\n');
    for i=1:N_clones
        fprintf('%d\n',P_conso(2*i));
    end
    fprintf('%d\n',P_conso(1));
    fprintf('Relapse clonal percentages from patient data:\n');
    disp(rel_clonal);
    fprintf('Relapse clonal percentages from the model:\n');
    disp(percent_BM_Relap);
    fileID = fopen(filename,'w');
        fprintf(fileID,'Leukemic parameters: (columns: a, p, d)\n');
        for j=1:N_clones
            fprintf(fileID,'%f ',parameters_l(:,j));
            fprintf(fileID,'\n');
        end
        fprintf(fileID,'\n\n\n');
        fprintf(fileID,'Population right after induction:\n');
        for j=1:length(P_induc)
            fprintf(fileID,'%f\n',P_induc(j));
        end
        fprintf(fileID,'\n\n\n');

        fprintf(fileID,'Population right before consolidation (biopsy):\n');
        for j=1:length(P_biopsy)
            fprintf(fileID,'%f\n',P_biopsy(j));
        end
        fprintf(fileID,'\n');
        fprintf(fileID,'Clonal percentages right before consolidation (biopsy):\n');
        for j=1:(N_clones+1)
            fprintf(fileID,'%f\n',percent_BM_biopsy(j));
        end
        fprintf(fileID,'\n\n\n');

        fprintf(fileID,'Population right after consolidation:\n');
        for j=1:length(P_conso)
            fprintf(fileID,'%f\n',P_conso(j));
        end
        fprintf(fileID,'\n');
        fprintf(fileID,'Clonal percentages right after consolidation:\n');
        for j=1:(N_clones+1)
            fprintf(fileID,'%f\n',percent_BM_Conso(j));
        end
        fprintf(fileID,'\n\n\n');
        fprintf(fileID,'Population at relapse:\n');
        for j=1:length(P_conso)
            fprintf(fileID,'%f\n',P_relap(j));
        end
        fprintf(fileID,'\n');
        fprintf(fileID,'Clonal percentages at relapse:\n');
        for j=1:(N_clones+1)
            fprintf(fileID,'%f\n',percent_BM_Relap(j));
        end
    fclose(fileID);
end
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
