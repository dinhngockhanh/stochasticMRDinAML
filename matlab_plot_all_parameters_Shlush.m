figure(1)
clf
for i_patient=1:7
%---Parameters for each patient
    if i_patient==1
        patient_ID          = 'Shlush2-1';
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
    elseif i_patient==2
        patient_ID          = 'Shlush2-2';
        index_parameter     = 1;
        N_clones            = 5;
        weight              = 86.9;
        p_cellularity       = 0.9;
        diag_PB_normal      = 3.00*(10^9);
        diag_PB_blast       = 2.79*(10^11);
        diag_BM_clonal      = [3.32     62.82   5.15    16.85   1.86    10.00];
        diag_PB_clonal      = [3.69     69.80   5.72    18.73   2.06];
        rel_BM_clonal       = [29.09    3.92    25.97   2.19    8.83    30.00];
        L_induction_1       = 7;
        k_induction_1       = [1*0.55    5*0.55];
        L_induction_2       = 0;
        k_induction_2       = [1*0.55    5*0.55];
        L_consolidation_1   = 10;
        k_consolidation_1   = [1*0.55    5*0.55];
        L_consolidation_2   = 10;
        k_consolidation_2   = [1*0.55    5*0.55];
        L_relapse           = 480;
    elseif i_patient==3
        patient_ID          = 'Shlush3-1';
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
    elseif i_patient==4
        patient_ID          = 'Shlush3-2';
        index_parameter     = 1;
        N_clones            = 4;
        weight              = 45.0;
        p_cellularity       = 0.9;
        diag_PB_normal      = 6.90*(10^10);
        diag_PB_blast       = 4.74*(10^10);
        diag_BM_clonal      = [20.26    19.61   9.05    1.08    50.00];
        diag_PB_clonal      = [40.52    39.23   18.10   2.16];
        rel_BM_clonal       = [7.53     22.31   0.93    4.23    65.00];
        L_induction_1       = 7;
        k_induction_1       = [1*0.9    4*0.9];
        L_induction_2       = 0;
        k_induction_2       = [1*0.9    4*0.9];
        L_consolidation_1   = 10;
        k_consolidation_1   = [1*0.9    4*0.9];
        L_consolidation_2   = 0;
        k_consolidation_2   = [1*0.9    4*0.9];
        L_relapse           = 180;
    elseif i_patient==5
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
    elseif i_patient==6
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
    elseif i_patient==7
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
    N_parameters        = 2*N_clones;
%---Parameters for the subplot
    % filename            = ['Patient modules and parameters/' patient_ID '.txt'];
    filename            = ['parameters_' patient_ID '.txt'];
    filename_1          = ['parameters_' patient_ID '_left.txt'];
    filename_2          = ['parameters_' patient_ID '_centered.txt'];
    filename_3          = ['parameters_' patient_ID '_right.txt'];
    if (exist(filename,'file')~=2)||(exist(filename_1,'file')~=2)||(exist(filename_2,'file')~=2)||(exist(filename_3,'file')~=2)
        continue;
    end

    fileID              = fopen(filename,'r');
    formatSpec          = '%f';
    PARAMETERS          = fscanf(fileID,formatSpec);
    PARAMETERS          = reshape(PARAMETERS,N_parameters,[])';
    fclose(fileID);
    N_parameter_sets    = size(PARAMETERS,1);

    min_a_l = min(min(PARAMETERS(:,1:N_clones)));
    max_a_l = max(max(PARAMETERS(:,1:N_clones)));
    min_p_l = min(min(PARAMETERS(:,N_clones+1:2*N_clones)));
    max_p_l = max(max(PARAMETERS(:,N_clones+1:2*N_clones)));

    min_a_l = 0.01*floor(100*min_a_l);
    max_a_l = 0.01*ceil(100*max_a_l);
    min_p_l = 0.01*floor(100*min_p_l);
    max_p_l = 0.01*ceil(100*max_p_l);

    fileID              = fopen(filename_1,'r');
    formatSpec          = '%f';
    parameter_left      = fscanf(fileID,formatSpec);
    fclose(fileID);

    fileID              = fopen(filename_2,'r');
    formatSpec          = '%f';
    parameter_centered  = fscanf(fileID,formatSpec);
    fclose(fileID);

    fileID              = fopen(filename_3,'r');
    formatSpec          = '%f';
    parameter_right     = fscanf(fileID,formatSpec);
    fclose(fileID);
%---Plot all parameter sets
    subplot(2,4,i_patient);
    for I_PARA=1:N_clones
        p = plot(PARAMETERS(1:1:end,I_PARA),PARAMETERS(1:1:end,I_PARA+N_clones),'o');hold on
        p.MarkerSize    = 5;
        p.LineWidth     = 1.5;
        if I_PARA==1
            p.MarkerFaceColor   = 'w';
            p.MarkerEdgeColor   = 'b';
        elseif I_PARA==2
            p.MarkerFaceColor = 'w';
            p.MarkerEdgeColor = 'm';
        elseif I_PARA==3
            p.MarkerFaceColor = 'w';
            p.MarkerEdgeColor = 'c';
        elseif I_PARA==4
            p.MarkerFaceColor = 'w';
            p.MarkerEdgeColor = 'k';
        elseif I_PARA==5
            p.MarkerFaceColor = 'w';
            p.MarkerEdgeColor = [0.45 0.45 0.45];
        end
        p.DisplayName   = ['Leukemic clone ' num2str(I_PARA)];
    end
%---Plot the three chosen parameter sets
    for j=1:N_clones
        p1  = plot(parameter_left(j),parameter_left(j+N_clones),'<');hold on
        p2  = plot(parameter_centered(j),parameter_centered(j+N_clones),'o');hold on
        p3  = plot(parameter_right(j),parameter_right(j+N_clones),'>');hold on
        p1.MarkerSize   = 15;
        p1.LineWidth    = 3;
        if j==1
            p1.MarkerFaceColor  = 'b';
        elseif j==2
            p1.MarkerFaceColor  = 'm';
        elseif j==3
            p1.MarkerFaceColor  = 'c';
        elseif j==4
            p1.MarkerFaceColor  = 'k';
        elseif j==5
            p1.MarkerFaceColor  = [0.45 0.45 0.45];
        end
        p1.MarkerEdgeColor      = 'k';
        p2.MarkerSize           = p1.MarkerSize;
        p3.MarkerSize           = p1.MarkerSize;
        p2.LineWidth            = p1.LineWidth;
        p3.LineWidth            = p1.LineWidth;
        p2.MarkerFaceColor      = p1.MarkerFaceColor;
        p3.MarkerFaceColor      = p1.MarkerFaceColor;
        p2.MarkerEdgeColor      = p1.MarkerEdgeColor;
        p3.MarkerEdgeColor      = p1.MarkerEdgeColor;
        set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        set(get(get(p2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        set(get(get(p3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
%---Configures for the plots
    xlim([min_a_l max_a_l]);
    ylim([min_p_l max_p_l]);
    t = title(['Patient ' num2str(i_patient) ' (ID: ' num2str(patient_ID) ')']);
    t.FontSize  = 15;
    l = xlabel('Renewal rate a_l');
    l.FontSize  = 15;
    if i_patient==1
        l = ylabel('Proliferation rate p_l');
        l.FontSize  = 15;
    end


    if i_patient==2
        lgd     = legend('show');
        lgd.FontSize    = 15;
    end
end

figure_name = 'plot_all_parameters_Shlush.png';
saveas(gcf,figure_name);
