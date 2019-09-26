figure(1)
clf
for i_patient=1:6
%---Parameters for each patient
    if i_patient==1
        patient_ID          = '400220';
        N_clones            = 2;
        weight              = 70.3;
        p_cellularity       = 0.9;
        diag_normal         = 2.94*(10^9);
        diag_blast          = 5.68*(10^10);
        diag_clonal         = [68   3   29]';
        rel_clonal          = [38   22  40]';
        L_induction         = 7;
        k_induction         = [0.8 3];
        L_consolidation     = 11;
        k_consolidation     = [1 2];
        L_relapse           = 252;
    elseif i_patient==2
        patient_ID          = '426980';
        N_clones            = 5;
        weight              = 87.5;
        p_cellularity       = 0.5;
        diag_normal         = 1.65*(10^9);
        diag_blast          = 1.32*(10^10);
        diag_clonal         = [3    3   33  22  3   36]';
        rel_clonal          = [0    3   0   0   9   88]';
        L_induction         = 7;
        k_induction         = [0.8  3.0];
        L_consolidation     = 8;
        k_consolidation     = [1.0  3.5];
        L_relapse           = 805;
    elseif i_patient==3
        patient_ID          = '452198';
        N_clones            = 4;
        weight              = 103.6;
        p_cellularity       = 0.9;
        diag_normal         = 3.63*(10^9);
        diag_blast          = 2.90*(10^10);
        diag_clonal         = [19   52  23  3   3]';
        rel_clonal          = [0    0   0   20  80]';
        L_induction         = 7;
        k_induction         = [0.6  3.0];
        L_consolidation     = 11;
        k_consolidation     = [1.0  4.0];
        L_relapse           = 505;
    elseif i_patient==4
        patient_ID          = '573988';
        N_clones            = 2;
        weight              = 80.7;
        p_cellularity       = 0.9;
        diag_normal         = 1.22*(10^10);
        diag_blast          = 7.60*(10^9);
        diag_clonal         = [72   3   25]';
        rel_clonal          = [15   39  46]';
        L_induction         = 7;
        k_induction         = [1.0  4.0];
        L_consolidation     = 5;
        k_consolidation     = [0.5  2.5];
        L_relapse           = 365;
    elseif i_patient==5
        patient_ID          = '758168';
        N_clones            = 3;
        weight              = 90.2;
        p_cellularity       = 0.9;
        diag_normal         = 3.26*(10^9);
        diag_blast          = 0.00*(10^1);
        diag_clonal         = [7    83  3   7]';
        rel_clonal          = [39   18  35  8]';
        L_induction         = 21;
        k_induction         = [0.5  1.0];
        L_consolidation     = 90;
        k_consolidation     = [0.75 1.0];
        L_relapse           = 961;
    elseif i_patient==6
        patient_ID          = '804168';
        N_clones            = 2;
        weight              = 79.3;
        p_cellularity       = 0.9;
        diag_normal         = 1.54*(10^11);
        diag_blast          = 2.29*(10^11);
        diag_clonal         = [83   3   14]';
        rel_clonal          = [36   45  19]';
        L_induction         = 7;
        k_induction         = [0.8  3.0];
        L_consolidation     = 11;
        k_consolidation     = [1.0  2.0];
        L_relapse           = 235;
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
    subplot(2,3,i_patient);
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

figure_name = 'plot_all_parameters.png';
saveas(gcf,figure_name);
