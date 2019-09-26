function matlab_plotting_deterministic
%-------------------------------------------------Options for the algorithm
    global a_c p_c d_c k C2 C1 m
%   Parameters of hematopoietic clone:
    a_c             = 0.87;
    p_c             = 0.45;
    d_c             = 2.3;
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

    N_patients      = 11;
    xlabel_ID       = cell(1,N_patients);

    vec_WBC = zeros(1,N_patients);

    figure(1);clf

    for i_patient=1:N_patients
        if i_patient==1
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
        elseif i_patient==2
            patient_ID          = '426980';
            index_parameter     = 1;
            N_clones            = 5;
            weight              = 87.5;
            p_cellularity       = 0.5;
            diag_normal         = 1.65*(10^9);
            diag_blast          = 1.32*(10^10);
            diag_clonal         = [3    3   33  22  3   36]';
            rel_clonal          = [0    3   0   0   9   88]';
            L_induction         = 7;
            k_induction         = [1*0.9 5*0.9];
            L_consolidation     = 8;
            k_consolidation     = [1*0.9 5*0.9];
            L_relapse           = 805;
        elseif i_patient==3
            patient_ID          = '452198';
            index_parameter     = 1;
            N_clones            = 4;
            weight              = 103.6;
            p_cellularity       = 0.9;
            diag_normal         = 3.63*(10^9);
            diag_blast          = 2.90*(10^10);
            diag_clonal         = [19   52  23  3   3]';
            rel_clonal          = [0    0   0   20  80]';
            L_induction         = 7;
            k_induction         = [1*1.0 5*1.0];
            L_consolidation     = 11;
            k_consolidation     = [1*1.0 5*1.0];
            L_relapse           = 505;
        elseif i_patient==4
            patient_ID          = '573988';
            index_parameter     = 1;
            N_clones            = 2;
            weight              = 80.7;
            p_cellularity       = 0.9;
            diag_normal         = 1.22*(10^10);
            diag_blast          = 7.60*(10^9);
            diag_clonal         = [72   3   25]';
            rel_clonal          = [15   39  46]';
            L_induction         = 7;
            k_induction         = [1*1.0 5*1.0];
            L_consolidation     = 5;
            k_consolidation     = [1*1.0 5*1.0];
            L_relapse           = 365;
        elseif i_patient==5
            patient_ID          = '758168';
            index_parameter     = 1;
            N_clones            = 3;
            weight              = 90.2;
            p_cellularity       = 0.9;
            diag_normal         = 3.26*(10^9);
            diag_blast          = 0.0;
            diag_clonal         = [7    83  3   7]';
            rel_clonal          = [39   18  35  8]';
            L_induction         = 21;
            k_induction         = [1*0.15 5*0.15];
            L_consolidation     = 90;
            k_consolidation     = [1*0.15 5*0.15];
            L_relapse           = 961;
        elseif i_patient==6
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
        elseif i_patient==7
            patient_ID          = 'Shlush2';
            index_parameter     = 1;
            N_clones            = 5;
            weight              = 86.9;
            p_cellularity       = 0.9;
            diag_PB_normal      = 3.00*(10^9);
            diag_PB_blast       = 2.79*(10^11);
            diag_BM_clonal      = [1.46     62.82   7.01    16.85   1.86    10.00];
            diag_PB_clonal      = [1.63     69.80   7.78    18.73   2.06];
            rel_BM_clonal       = [20.26    3.92    34.80   2.19    8.83    30.00];
            L_induction_1       = 7;
            k_induction_1       = [1*0.55    5*0.55];
            L_induction_2       = 0;
            k_induction_2       = [1*0.55    5*0.55];
            L_consolidation_1   = 10;
            k_consolidation_1   = [1*0.55    5*0.55];
            L_consolidation_2   = 10;
            k_consolidation_2   = [1*0.55    5*0.55];
            L_relapse           = 480;
        elseif i_patient==8
            patient_ID          = 'Shlush3';
            index_parameter     = 1;
            N_clones            = 4;
            weight              = 45.0;
            p_cellularity       = 0.9;
            diag_PB_normal      = 6.90*(10^10);
            diag_PB_blast       = 4.74*(10^10);
            diag_BM_clonal      = [19.18    20.69   9.05    1.08    50.00];
            diag_PB_clonal      = [38.36    41.39   18.10   2.16];
            rel_BM_clonal       = [3.30     26.54   0.93    4.23    65.00];
            L_induction_1       = 7;
            k_induction_1       = [1*0.9    4*0.9];
            L_induction_2       = 0;
            k_induction_2       = [1*0.9    4*0.9];
            L_consolidation_1   = 10;
            k_consolidation_1   = [1*0.9    4*0.9];
            L_consolidation_2   = 0;
            k_consolidation_2   = [1*0.9    4*0.9];
            L_relapse           = 180;
        elseif i_patient==9
            patient_ID          = 'Shlush6';
            index_parameter     = 1;
            N_clones            = 5;
            weight              = 76.0;
            p_cellularity       = 0.9;
            diag_PB_normal      = 3.80*(10^10);
            diag_PB_blast       = 6.48*(10^10);
            diag_BM_clonal      = [1.13     4.14    34.09   0.93    9.71    50.00];
            diag_PB_clonal      = [2.27     8.29    68.17   1.85    19.42];
            rel_BM_clonal       = [0.00     33.51   0.00    9.49    0.00    57.00];
            L_induction_1       = 7;
            k_induction_1       = [1*0.7    4*0.7];
            L_induction_2       = 0;
            k_induction_2       = [1*0.7    4*0.7];
            L_consolidation_1   = 10;
            k_consolidation_1   = [1*0.7    4*0.7];
            L_consolidation_2   = 10;
            k_consolidation_2   = [1*0.7    4*0.7];
            L_relapse           = 480;
        elseif i_patient==10
            patient_ID          = 'Shlush8';
            index_parameter     = 1;
            N_clones            = 2;
            weight              = 90.0;
            p_cellularity       = 0.9;
            diag_PB_normal      = 2.00*(10^9);
            diag_PB_blast       = 1.69*(10^11);
            diag_BM_clonal      = [78.47    8.53    13.00];
            diag_PB_clonal      = [90.20    9.80];
            rel_BM_clonal       = [82.26    12.74   5.00];
            L_induction_1       = 7;
            k_induction_1       = [1*0.7    4*0.7];
            L_induction_2       = 7;
            k_induction_2       = [1*0.7    4*0.7];
            L_consolidation_1   = 10;
            k_consolidation_1   = [1*0.7    4*0.7];
            L_consolidation_2   = 10;
            k_consolidation_2   = [1*0.7    4*0.7];
            L_relapse           = 240;
        elseif i_patient==11
            patient_ID          = 'Shlush15';
            index_parameter     = 1;
            N_clones            = 4;
            weight              = 85.0;
            p_cellularity       = 0.9;
            diag_PB_normal      = 3.04*(10^10);
            diag_PB_blast       = 1.02*(10^10);
            diag_BM_clonal      = [70.65    0.03    16.98   1.34    11.00];
            diag_PB_clonal      = [79.38    0.03    19.07   1.52];
            rel_BM_clonal       = [0.00     47.47   3.52    10.01   39.00];
            L_induction_1       = 7;
            k_induction_1       = [1*0.7    4*0.7];
            L_induction_2       = 0;
            k_induction_2       = [1*0.7    4*0.7];
            L_consolidation_1   = 10;
            k_consolidation_1   = [1*0.7    4*0.7];
            L_consolidation_2   = 10;
            k_consolidation_2   = [1*0.7    4*0.7];
            L_relapse           = 300;
        end
        xlabel_ID{i_patient}    = patient_ID;
        if (i_patient==7)||(i_patient==8)
            patient_ID  = [patient_ID '-1']
        end
        if i_patient<=6
            [WBC_induc,WBC_conso] = WBC_nadirs_Ding(patient_ID,N_clones,weight,...
                                                    p_cellularity,diag_clonal,...
                                                    diag_normal,diag_blast,...
                                                    L_induction,k_induction,...
                                                    L_consolidation,k_consolidation,...
                                                    L_relapse);
        else
            [WBC_induc,WBC_conso] = WBC_nadirs_Shlush(patient_ID,N_clones,weight,...
                                                      p_cellularity,diag_BM_clonal,...
                                                      diag_PB_normal,diag_PB_blast,diag_PB_clonal,...
                                                      L_induction_1,k_induction_1,...
                                                      L_induction_2,k_induction_2,...
                                                      L_consolidation_1,k_consolidation_1,...
                                                      L_consolidation_2,k_consolidation_2,...
                                                      L_consolidation,k_consolidation,...
                                                      L_relapse)
        end
        vec_WBC(i_patient)  = min(WBC_induc,WBC_conso);
    end

    figure(1)
    if N_patients<=6
        p = scatter([1:N_patients],vec_WBC,200);
        p.Marker            = 'o';
        p.MarkerEdgeColor   = 'b';
        p.MarkerFaceColor   = 'b';
    else
        p = scatter([1:6],vec_WBC(1:6),200);hold on
        p.Marker            = 'o';
        p.MarkerEdgeColor   = 'b';
        p.MarkerFaceColor   = 'b';
        p = scatter([7:N_patients],vec_WBC(7:end),200);hold on
        p.Marker            = 'o';
        p.MarkerEdgeColor   = 'm';
        p.MarkerFaceColor   = 'm';
    end

    plot([0,N_patients+1],[0.1,0.1],'r--','LineWidth',3)
    plot([0,N_patients+1],[1,1],'r--','LineWidth',3)

    xlabel('Patient');
    ylabel('WBC nadir (10^9 cells/L)')

    set(gca,'fontsize',20);
    xticklabels(xlabel_ID);xtickangle(45);

    set(gca,'xtick',[1:N_patients]);

    figure_name = ['plot_WBC_nadir.png'];
    saveas(gcf,figure_name);

end

function [WBC_induc,WBC_conso] = WBC_nadirs_Shlush(patient_ID,N_clones,weight,...
                                                   p_cellularity,diag_BM_clonal,...
                                                   diag_PB_normal,diag_PB_blast,diag_PB_clonal,...
                                                   L_induction_1,k_induction_1,...
                                                   L_induction_2,k_induction_2,...
                                                   L_consolidation_1,k_consolidation_1,...
                                                   L_consolidation_2,k_consolidation_2,...
                                                   L_consolidation,k_consolidation,...
                                                   L_relapse)
    global a_c p_c d_c k C2 C1 m
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
%--------------------------------------------Set up parameters at diagnosis
    parameters_l= [para_a_l;  para_p_l;  para_d_l];
%  	Population vector at diagnosis
    P           = zeros(1,2*N_clones+2);
    P(1)        = diag_BM_sum*diag_BM_clonal(end)/100;
    P(2)        = diag_PB_normal;
    for j=1:N_clones
%       Compute number of mitotic cells in BM
        P(2*j+1)= diag_BM_sum*diag_BM_clonal(j)/100;
%       Compute number of mature cells in blood
        P(2*j+2)= diag_PB_blast*diag_PB_clonal(j)/100;
    end
%-------------------------------------From diagnosis to end of chemotherapy
%   Run the ODEs for induction 1
    T_temp  = 0;
    ODE     = @(t,P)ODE_chemo(t,P,parameters_l,m,k,C1,C2,...
                                  k_induction_1,c_1_bar,a_c,p_c,d_c);
    [T_result,P_result]=ode23(ODE,[T_temp (T_temp+L_induction_1)],P);
    T_temp  = T_temp+L_induction_1;
    T_entire= T_result;
    P_entire= P_result;
    P_temp  = P_result(end,:);
    if L_induction_2>0
%       Run the ODEs for the gap between induction 1 and induction 2
        ODE     = @(t,P)ODE_normal(t,P,parameters_l,m,k,C1,C2,...
                                   c_1_bar,a_c,p_c,d_c);
        [T_result,P_result]=ode23(ODE,[T_temp (T_temp+3)],P_temp);
%       Save the result from the ODEs
        T_entire= [T_entire' T_result']';
        P_entire= [P_entire' P_result']';
        T_temp  = T_temp+3;
        P_temp  = P_result(end,:);
%       Run the ODEs for induction 2
        ODE     = @(t,P)ODE_chemo(t,P,parameters_l,m,k,C1,C2,...
                                      k_induction_2,c_1_bar,a_c,p_c,d_c);
        [T_result,P_result]=ode23(ODE,[T_temp (T_temp+L_induction_2)],P_temp);
        T_temp  = T_temp+L_induction_2;
%       Save the result from the ODEs
        T_entire= [T_entire' T_result']';
        P_entire= [P_entire' P_result']';
    end
    P_induc = P_result(end,:);
    P_temp  = P_result(end,:);
%   Run the ODEs for the gap between induction and consolidation
    ODE     = @(t,P)ODE_normal(t,P,parameters_l,m,k,C1,C2,...
                               c_1_bar,a_c,p_c,d_c);
    [T_result,P_result]=ode23(ODE,[T_temp (T_temp+30)],P_temp);
    T_temp  = T_temp+30;
    T_entire= [T_entire' T_result']';
    P_entire= [P_entire' P_result']';
    P_biopsy= P_result(end,:);
    P_temp  = P_result(end,:);
%   Calculate clonal and normal percentages in BM
    percent_BM_biopsy    = zeros(N_clones+1,1);
    Sum_BM  = sum(P_biopsy(3:2:end))+P_biopsy(1);
    for j=1:N_clones
        percent_BM_biopsy(j)=100*P_biopsy(2*j+1)/Sum_BM;
    end
    percent_BM_biopsy(N_clones+1)=100*P_biopsy(1)/Sum_BM;
%   Run the ODEs for consolidation 1
    ODE     = @(t,P)ODE_chemo(t,P,parameters_l,m,k,C1,C2,...
                              k_consolidation_1,c_1_bar,a_c,p_c,d_c);
    [T_result,P_result]=ode23(ODE,[T_temp (T_temp+L_consolidation_1)],P_temp);
    T_temp  = T_temp+L_consolidation_1;
    T_entire= [T_entire' T_result']';
    P_entire= [P_entire' P_result']';
    P_conso = P_result(end,:);
    P_temp  = P_result(end,:);
    if L_consolidation_2>0
%       Run the ODEs for the gap between consolidation 1 and consolidation 2
        ODE     = @(t,P)ODE_normal(t,P,parameters_l,m,k,C1,C2,...
                                   c_1_bar,a_c,p_c,d_c);
        [T_result,P_result]=ode23(ODE,[T_temp (T_temp+10)],P_temp);
        T_temp  = T_temp+10;
        P_temp  = P_result(end,:);
%       Save the result from the ODEs
        T_entire= [T_entire' T_result']';
        P_entire= [P_entire' P_result']';
%       Run the ODEs for consolidation 2
        ODE     = @(t,P)ODE_chemo(t,P,parameters_l,m,k,C1,C2,...
                                  k_consolidation_2,c_1_bar,a_c,p_c,d_c);
        [T_result,P_result]=ode23(ODE,[T_temp (T_temp+L_consolidation_2)],P_temp);
        T_temp  = T_temp+L_consolidation_2;
        P_conso = P_result(end,:);
        P_temp  = P_result(end,:);
%       Save the result from the ODEs
        T_entire= [T_entire' T_result']';
        P_entire= [P_entire' P_result']';
    end
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
%   Calculate clonal and normal percentages in PB
    percent_PB_Relap    = zeros(N_clones+1,1);
    Sum_PB              = sum(P_relap(2:2:end));
    for j=1:N_clones
        percent_PB_Relap(j)=100*P_relap(2*j+2)/Sum_PB;
    end
    percent_PB_Relap(N_clones+1)=100*P_relap(2)/Sum_PB;
%--------------------------------------------------------Output the outcome
%   WBC at end of induction
    WBC_induc   = sum(P_induc(2:2:end))/5/10^9;
    WBC_conso   = sum(P_conso(2:2:end))/5/10^9;
end

function [WBC_induc,WBC_conso] = WBC_nadirs_Ding(patient_ID,N_clones,weight,...
                                                 p_cellularity,diag_clonal,...
                                                 diag_normal,diag_blast,...
                                                 L_induction,k_induction,...
                                                 L_consolidation,k_consolidation,...
                                                 L_relapse)
    global a_c p_c d_c k C2 C1 m
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
%--------------------------------------------------------Output the outcome
%   WBC at end of induction
    WBC_induc   = sum(P_induc(2:2:end))/5/10^9;
    WBC_conso   = sum(P_conso(2:2:end))/5/10^9;
end
%==========================================================================
function dP = ODE_normal(t,P,parameters_l,m,k,C1,C2,c_1_bar,a_c,p_c,d_c)
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
function dP = ODE_chemo(t,P,parameters_l,m,k,C1,C2,...
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
