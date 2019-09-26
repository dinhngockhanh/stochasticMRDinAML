%-----------------------------------------------Parameters for each patient
% patient_ID          = '400220';
% N_clones            = 2;
% weight              = 70.3;
% p_cellularity       = 0.9;
% diag_normal         = 2.94*(10^9);
% diag_blast          = 5.68*(10^10);
% diag_clonal         = [68   3   29]';
% rel_clonal          = [38   22  40]';
% L_induction         = 7;
% k_induction         = [0.8 1];
% L_consolidation     = 11;
% k_consolidation     = [0.3 0.2];
% L_relapse           = 252;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% patient_ID          = '426980';
% N_clones            = 5;
% weight              = 87.5;
% p_cellularity       = 0.5;
% diag_normal         = 1.65*(10^9);
% diag_blast          = 1.32*(10^10);
% diag_clonal         = [3    3   33  22  3   36]';
% rel_clonal          = [0    3   0   0   9   88]';
% L_induction         = 7;
% k_induction         = [0.8  3.0];
% L_consolidation     = 8;
% k_consolidation     = [1.0  3.5];
% L_relapse           = 805;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% patient_ID          = '573988';
% N_clones            = 2;
% weight              = 80.7;
% p_cellularity       = 0.9;
% diag_normal         = 1.22*(10^10);
% diag_blast          = 7.60*(10^9);
% diag_clonal         = [72   3   25]';
% rel_clonal          = [15   39  46]';
% L_induction         = 7;
% k_induction         = [1.0  4.0];
% L_consolidation     = 5;
% k_consolidation     = [0.5  2.5];
% L_relapse           = 365;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% patient_ID          = '758168';
% N_clones            = 3;
% weight              = 90.2;
% p_cellularity       = 0.9;
% diag_normal         = 3.26*(10^9);
% diag_blast          = 0.00*(10^1);
% diag_clonal         = [7    83  3   7]';
% rel_clonal          = [39   18  35  8]';
% L_induction         = 21;
% k_induction         = [0.5  1.0];
% L_consolidation     = 90;
% k_consolidation     = [0.75 1.0];
% L_relapse           = 961;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% patient_ID          = '804168';
% N_clones            = 2;
% weight              = 79.3;
% p_cellularity       = 0.9;
% diag_normal         = 1.54*(10^11);
% diag_blast          = 2.29*(10^11);
% diag_clonal         = [83   3   14]';
% rel_clonal          = [36   45  19]';
% L_induction         = 7;
% k_induction         = [0.8  3.0];
% L_consolidation     = 11;
% k_consolidation     = [1.0  2.0];
% L_relapse           = 235;





%=======================================================================
min_a_l	= 0.7;
max_a_l	= 1.0;
min_p_l	= 0.3;
max_p_l	= 0.7;

center_a_l  = 0.85;
center_p_l  = 0.5;

left_a_l    = 0.7;
left_p_l    = 0.7;

right_a_l   = 1;
right_p_l   = 0.3;

% filename            = ['Patient modules and parameters/' patient_ID '.txt'];
filename            = ['parameters_' patient_ID '.txt'];
N_parameters        = 2*N_clones;

fileID              = fopen(filename,'r');
formatSpec          = '%f';
PARAMETERS          = fscanf(fileID,formatSpec);
PARAMETERS          = reshape(PARAMETERS,N_parameters,[])';
fclose(fileID);
N_parameter_sets    = size(PARAMETERS,1);

%============Find the "left", "center" and "right" of the parameter sets
min_a_l_real    = min(min(PARAMETERS(:,1:N_clones)));
max_a_l_real    = max(max(PARAMETERS(:,1:N_clones)));
min_p_l_real    = min(min(PARAMETERS(:,N_clones+1:2*N_clones)));
max_p_l_real    = max(max(PARAMETERS(:,N_clones+1:2*N_clones)));

left_a_l        = min_a_l_real;
left_p_l        = max_p_l_real;

right_a_l       = max_a_l_real;
right_p_l       = min_p_l_real;

center_a_l      = 0.5*(left_a_l+right_a_l);
center_p_l      = 0.5*(left_p_l+right_p_l);
%=========================Find the parameter set closest to the "center"
vec_error_from_center   = zeros(1,N_parameter_sets);
for i=1:N_parameter_sets
    parameter       = PARAMETERS(i,:);
    error           = 0;
    for j=1:N_clones
        error       = error+abs(parameter(j)-center_a_l)+...
                            abs(parameter(j+N_clones)-center_p_l);
        vec_error_from_center(i)    = error;
    end
end
[~,pos]             = min(vec_error_from_center);
fprintf('Centered at position %d\n',pos);
parameter_centered  = PARAMETERS(pos,:);
fileID  = fopen(['parameters_' patient_ID '_centered.txt'],'w');
for i=1:N_parameters
    fprintf(fileID,'%1.16f ',parameter_centered(i));
end
fprintf(fileID,'\n');
fclose(fileID);
%===========================Find the parameter set closest to the "left"
vec_error_from_left = zeros(1,N_parameter_sets);
for i=1:N_parameter_sets
    parameter       = PARAMETERS(i,:);
    error           = 0;
    for j=1:N_clones
        error       = error+abs(parameter(j)-left_a_l)+...
                            abs(parameter(j+N_clones)-left_p_l);
        vec_error_from_left(i)  = error;
    end
end
[~,pos]             = min(vec_error_from_left);
fprintf('Left at position %d\n',pos);
parameter_left      = PARAMETERS(pos,:);
fileID  = fopen(['parameters_' patient_ID '_left.txt'],'w');
for i=1:N_parameters
    fprintf(fileID,'%1.16f ',parameter_left(i));
end
fprintf(fileID,'\n');
fclose(fileID);
%==========================Find the parameter set closest to the "right"
vec_error_from_right = zeros(1,N_parameter_sets);
for i=1:N_parameter_sets
    parameter       = PARAMETERS(i,:);
    error           = 0;
    for j=1:N_clones
        error       = error+abs(parameter(j)-right_a_l)+...
                            abs(parameter(j+N_clones)-right_p_l);
        vec_error_from_right(i)  = error;
    end
end
[~,pos]             = min(vec_error_from_right);
fprintf('Right at position %d\n',pos);
parameter_right      = PARAMETERS(pos,:);
fileID  = fopen(['parameters_' patient_ID '_right.txt'],'w');
for i=1:N_parameters
    fprintf(fileID,'%1.16f ',parameter_right(i));
end
fprintf(fileID,'\n');
fclose(fileID);
%================================================Plot the parameter sets
figure(1)
clf
for i_fig=1:3
    subplot(3,1,i_fig);
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
            p.MarkerEdgeColor = 'k';
        elseif I_PARA==4
            p.MarkerFaceColor = 'w';
            p.MarkerEdgeColor = 'c';
        elseif I_PARA==5
            p.MarkerFaceColor = 'w';
            p.MarkerEdgeColor = [0.45 0.45 0.45];
        end
        p.DisplayName   = ['Cancer clone ' num2str(I_PARA)];
    end
    if i_fig==1
        lgd     = legend('show');
        lgd.FontSize    = 15;
    end
    for j=1:N_clones
        if i_fig==1
            p1  = plot(parameter_left(j),parameter_left(j+N_clones),'o');hold on
        elseif i_fig==2
            p1  = plot(parameter_centered(j),parameter_centered(j+N_clones),'o');hold on
        else
            p1  = plot(parameter_right(j),parameter_right(j+N_clones),'o');hold on
        end
        p1.MarkerSize    = 15;
        p1.LineWidth     = 3;
        if j==1
            p1.MarkerFaceColor   = 'b';
        elseif j==2
            p1.MarkerFaceColor   = 'm';
        elseif j==3
            p1.MarkerFaceColor   = 'k';
        elseif j==4
            p1.MarkerFaceColor   = 'c';
        elseif j==5
            p1.MarkerFaceColor   = [0.45 0.45 0.45];
        end
        p1.MarkerEdgeColor      = 'k';
    end

end

% t = text(PARAMETERS(1,1),PARAMETERS(1,5),num2str(1));
% t.FontSize  = 20;
% t.Color     = 'r';
% t = text(PARAMETERS(1,2),PARAMETERS(1,6),num2str(2));

for i_fig=1:3
    subplot(3,1,i_fig);
    ylim([min_p_l max_p_l]);
    xlim([min_a_l max_a_l]);
    t = title(['Parameter set ' num2str(i_fig) ' for patient ' patient_ID]);
    t.FontSize  = 20;
    set(gca,'FontSize',20);
end
subplot(3,1,3);
    l   = xlabel('Renewal rate a_l');
    l.FontSize  = 20;
subplot(3,1,2);
    l   = ylabel('Proliferation rate p_l');
    l.FontSize  = 20;

% return

% %-------------------------------Plots the ODE result for each parameter set
% for INDEX = 1:100
%     figure(2)
%     clf
% %   Parameters of hematopoietic clone:
%     a_c             = 0.87;
%     p_c             = 0.45;
%     d_c             = 3.2;
% %   The specific leukemic clones' parameters:
%     para_a_l       = PARAMETERS(INDEX,1:N_clones);
%     para_p_l       = PARAMETERS(INDEX,1+N_clones:2*N_clones);
%     para_d_l       = 2.5*ones(1,N_clones);
% %   Equilibrium hematopoietic mitotic cell numbers:
%     c_1_bar         = 2*(10^9)*weight;
% %   Total number of cells in bone marrow at primary/diagnosis:
%     full_BM_sum     = 4.6*(10^11);
%     diag_BM_sum     = p_cellularity*full_BM_sum;
% %   Constant for negative feedback of hematopoietic clone, as in
% %   feedback = 1/(1+k*c2) where c2 is number of post-mitotic hematopoietic
% %   cells number:
%     k               = 10^-12;
% %   Constants for negative feedback of overcrowding bone marrow, as in
% %   feedback = C1*max(0,x-C2*c_1_d_bar) where x is total number of cells in
% %   bone marrow (including mitotic hematopoietic and mitotic and 1st stage
% %   post-mitotic cancer cells):
%     C2              = 1;
%     C1              = 10^-12;
% %   Migration rate of post-mitotic cancer cells to peripheral blood:
%     m               = 10;
% %-----------------------------------------------Parameters of the algorithm
% %   Number of cancer clones there are:
%     N_cancer        = length(para_a_l);
% %--------------------------------------------Set up parameters at diagnosis
%     parameters_l= [para_a_l;  para_p_l;  para_d_l];
% %  	Population vector at diagnosis
%     P           = zeros(1,3*N_cancer+2);
%     P(1)        = diag_BM_sum*diag_clonal(end)/100;
%     P(2)        = diag_normal;
%     for j=1:N_cancer
% %     	Compute number of mitotic cells in BM
%         P(3*j)  = diag_BM_sum*diag_clonal(j)/100;
% %     	Compute number of mature cells in blood
%         P(3*j+2)= diag_blast*diag_clonal(j)/100;
%     end
% %-------------------------------------From diagnosis to end of chemotherapy
% %   Run the ODEs for induction
%     T_temp  = 0;
%     ODE     = @(t,P)ODE_chemo(t,P,parameters_l,m,k,C1,C2,...
%                                   k_induction,c_1_bar,a_c,p_c,d_c);
%     [T_result,P_result]=ode23(ODE,[T_temp (T_temp+L_induction)],P);
%     T_temp  = T_temp+L_induction;
%     P_induc = P_result(end,:);
%     P_temp  = P_result(end,:);
% %   Save the result from the ODEs
%     T_entire= T_result;
%     P_entire= P_result;
% %   Run the ODEs for the gap between induction and consolidation
%     ODE     = @(t,P)ODE_normal(t,P,parameters_l,m,k,C1,C2,...
%                                c_1_bar,a_c,p_c,d_c);
%     [T_result,P_result]=ode23(ODE,[T_temp (T_temp+30)],P_temp);
%     T_temp  = T_temp+30;
%     P_biopsy= P_result(end,:);
%     P_temp  = P_result(end,:);
% %   Save the result from the ODEs
%     T_entire= [T_entire' T_result']';
%     P_entire= [P_entire' P_result']';
% %   Calculate clonal and normal percentages in BM
%     percent_BM_Biopsy    = zeros(N_cancer+1,1);
%     Sum_BM  = sum(P_biopsy(3:3:end))+sum(P_biopsy(4:3:end))+P_biopsy(1);
%     for j=1:N_cancer
%         percent_BM_Biopsy(j)=100*(P_biopsy(3*j)+P_biopsy(3*j+1))/Sum_BM;
%     end
%     percent_BM_Biopsy(N_cancer+1)=100*P_biopsy(1)/Sum_BM;
% %   Run the ODEs for consolidation
%     ODE     = @(t,P)ODE_chemo(t,P,parameters_l,m,k,C1,C2,...
%                               k_consolidation,c_1_bar,a_c,p_c,d_c);
%     [T_result,P_result]=ode23(ODE,[T_temp (T_temp+L_consolidation)],P_temp);
%     T_temp  = T_temp+L_consolidation;
%     P_conso = P_result(end,:);
%     P_temp  = P_result(end,:);
% %   Save the result from the ODEs
%     T_entire= [T_entire' T_result']';
%     P_entire= [P_entire' P_result']';
% %   Calculate clonal and normal percentages in BM
%     percent_BM_Conso     = zeros(N_cancer+1,1);
%     Sum_BM  = sum(P_conso(3:3:end))+sum(P_conso(4:3:end))+P_conso(1);
%     for j=1:N_cancer
%         percent_BM_Conso(j)=100*(P_conso(3*j)+P_conso(3*j+1))/Sum_BM;
%     end
%     percent_BM_Conso(N_cancer+1)=100*P_conso(1)/Sum_BM;
% %---------------------------------------From end of chemotherapy to relapse
% %   Run the ODEs
%     ODE     = @(t,P)ODE_normal(t,P,parameters_l,m,k,C1,C2,...
%                                c_1_bar,a_c,p_c,d_c);
%     [T_result,P_result]=ode23(ODE,[T_temp L_relapse],P_temp);
%     P_relap = P_result(end,:);
% %   Save the result from the ODEs
%     T_entire= [T_entire' T_result']';
%     P_entire= [P_entire' P_result']';
% %   Calculate clonal and normal percentages in BM
%     percent_BM_Relap     = zeros(N_cancer+1,1);
%     Sum_BM  = sum(P_relap(3:3:end))+sum(P_relap(4:3:end))+P_relap(1);
%     for j=1:N_cancer
%         percent_BM_Relap(j)=100*(P_relap(3*j)+P_relap(3*j+1))/Sum_BM;
%     end
%     percent_BM_Relap(N_cancer+1)=100*P_relap(1)/Sum_BM;
% %-------------------------------------------------Plot the clonal evolution
% %   Some parameters for the plots
%     plot_1_max  = 1.2*max(max(P_entire(:,1)),max(max(P_entire(:,3:3:end))));
%     plot_2_max  = 1.2*max(max(P_entire(:,2)),max(max(P_entire(:,5:3:end))));
%     plot_1_min  = 1.2*min(min(P_entire(:,1)),min(min(P_entire(:,3:3:end))));
%     plot_2_min  = 1.2*min(min(P_entire(:,2)),min(min(P_entire(:,5:3:end))));
% %   Plot the clonal evolution in semilogy
%     subplot(2,20,[2:19]);     p = semilogy(T_entire,P_entire(:,1)); hold on
%                             p.LineWidth     = 2;
%                             p.Color         = 'r';
%     for j=1:N_cancer
%         p = semilogy(T_entire,P_entire(:,3*j)); hold on
%             p.LineWidth     = 2;
%             if j==1     p.Color     = 'b';
%             elseif j==2 p.Color     = 'm';
%             elseif j==3 p.Color     = 'k';
%             elseif j==4 p.Color     = 'c';
%             elseif j==5 p.Color     = [0.45 0.45 0.45];
%             end
%     end
%     p = rectangle('Position',[0 1 L_induction 10^ceil(log10(plot_1_max))],'FaceColor','g');hold on
%         uistack(p,'bottom');
%         p = rectangle('Position',[(L_induction+30) 1 L_consolidation 10^ceil(log10(plot_1_max))],'FaceColor','g');hold on
%         uistack(p,'bottom');
%         xlim([0 L_relapse]);    ylim([max(1,10^floor(log10(plot_1_min))) 10^ceil(log10(plot_1_max))]);
%         title(['Evolution of clonal populations in bone marrow (Parameter set ' num2str(INDEX) ')' ]);
%         vecYTick = [floor(log10(plot_1_min)):1:ceil(log10(plot_1_max))];
%         for i=1:length(vecYTick)
%             vecYTick(i) = 10^vecYTick(i);
%         end
%         ax = gca;
%         ax.YTick = vecYTick;
% %   Find the bone marrow clonal percentages and cellularity over time
%     for i=1:length(T_entire)
%         BM_cellular(i)  = P_entire(i,1)+sum(P_entire(i,3:3:end))+sum(P_entire(i,4:3:end));
%         for j=1:N_cancer
%             BM_clonal(i,j) = 100*(P_entire(i,3*j)+P_entire(i,3*j+1))/BM_cellular(i);
%         end
%         BM_clonal(i,N_cancer+1)  = 100*P_entire(i,1)/BM_cellular(i);
%         BM_cellular(i)  = 100*BM_cellular(i)/full_BM_sum;
%     end
%     max_BM_cellular = 20*ceil(max(BM_cellular)/20);
% %   Plot the clonal percentages over time
%     subplot(2,20,[22:39]);     p = area(BM_clonal);
%                         set(p,'XData',T_entire);
%                         p(N_cancer+1).FaceColor     = 'r';
%                         for j=1:N_cancer
%                             if j==1     p(j).FaceColor  = 'b';
%                             elseif j==2 p(j).FaceColor  = 'm';
%                             elseif j==3 p(j).FaceColor  = 'k';
%                             elseif j==4 p(j).FaceColor  = 'c';
%                             elseif j==5 p(j).FaceColor  = [0.45 0.45 0.45];
%                             end
%                         end
%                         set(gca,'ytick',[])
%                         xlim([0 L_relapse]);    ylim([0 100]);
%                         title('Evolution of clonal percentages');
% %   Plot the diagnosis and relapse clonal percentages
%     subplot(2,20,21); b = bar([diag_clonal';zeros(1,N_cancer+1)],'stacked');
%                        xlim([0.75 1.25]);
%                        b(N_cancer+1).FaceColor = 'r';
%                        for j=1:N_cancer
%                            if j==1      b(j).FaceColor  = 'b';
%                             elseif j==2 b(j).FaceColor  = 'm';
%                             elseif j==3 b(j).FaceColor  = 'k';
%                             elseif j==4 b(j).FaceColor  = 'c';
%                             elseif j==5 b(j).FaceColor  = [0.45 0.45 0.45];
%                             end
%                        end
%                        set(gca,'xtick',[]);
%     subplot(2,20,40); b = bar([rel_clonal';zeros(1,N_cancer+1)],'stacked');
%                        xlim([0.75 1.25]);
%                        b(N_cancer+1).FaceColor = 'r';
%                        for j=1:N_cancer
%                            if j==1      b(j).FaceColor  = 'b';
%                             elseif j==2 b(j).FaceColor  = 'm';
%                             elseif j==3 b(j).FaceColor  = 'k';
%                             elseif j==4 b(j).FaceColor  = 'c';
%                             elseif j==5 b(j).FaceColor  = [0.45 0.45 0.45];
%                             end
%                        end
%                        ax = gca;
%                        ax.YAxisLocation = 'right';
%                        set(gca,'xtick',[]);
% % %--------------------------------------------------------Output the outcome
% % %   Output results on screen
% %     fprintf('Clonal percentages right before consolidation (biopsy):\n');
% %     disp(percent_BM_Biopsy);
% %     fprintf('BM populations right after consolidation:\n');
% %     for i=1:N_cancer
% %         fprintf('%d\n',P_conso(3*i));
% %     end
% %     fprintf('%d\n',P_conso(1));
% %     fprintf('Relapse clonal percentages from patient data:\n');
% %     disp(rel_clonal);
% %     fprintf('Relapse clonal percentages from the model:\n');
% %     disp(percent_BM_Relap);
% %     fileID = fopen(filename,'w');
% %         fprintf(fileID,'Leukemic parameters: (columns: a, p, d)\n');
% %         for j=1:N_cancer
% %             fprintf(fileID,'%f ',parameters_l(:,j));
% %             fprintf(fileID,'\n');
% %         end
% %         fprintf(fileID,'\n\n\n');
% %         fprintf(fileID,'Population right after induction:\n');
% %         for j=1:length(P_induc)
% %             fprintf(fileID,'%f\n',P_induc(j));
% %         end
% %         fprintf(fileID,'\n\n\n');
% %
% %         fprintf(fileID,'Population right before consolidation (biopsy):\n');
% %         for j=1:length(P_biopsy)
% %             fprintf(fileID,'%f\n',P_biopsy(j));
% %         end
% %         fprintf(fileID,'\n');
% %         fprintf(fileID,'Clonal percentages right before consolidation (biopsy):\n');
% %         for j=1:(N_cancer+1)
% %             fprintf(fileID,'%f\n',percent_BM_Biopsy(j));
% %         end
% %         fprintf(fileID,'\n\n\n');
% %
% %         fprintf(fileID,'Population right after consolidation:\n');
% %         for j=1:length(P_conso)
% %             fprintf(fileID,'%f\n',P_conso(j));
% %         end
% %         fprintf(fileID,'\n');
% %         fprintf(fileID,'Clonal percentages right after consolidation:\n');
% %         for j=1:(N_cancer+1)
% %             fprintf(fileID,'%f\n',percent_BM_Conso(j));
% %         end
% %         fprintf(fileID,'\n\n\n');
% %         fprintf(fileID,'Population at relapse:\n');
% %         for j=1:length(P_conso)
% %             fprintf(fileID,'%f\n',P_relap(j));
% %         end
% %         fprintf(fileID,'\n');
% %         fprintf(fileID,'Clonal percentages at relapse:\n');
% %         for j=1:(N_cancer+1)
% %             fprintf(fileID,'%f\n',percent_BM_Relap(j));
% %         end
% %     fclose(fileID);
%
%     pause(2);
% end
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
    Sum         = sum(P(5:3:3*N_clones+2))+P(2);
    s           = 1/(1+k*Sum);
    x           = sum(P(3:3:3*N_clones+2))+sum(P(4:3:3*N_clones+2))+P(1);
    d           = C1*max(0,x-C2*c_1_bar);
%-------------------------------------------------------Deterministic model
%   Hematopoietic clone: (if mitotic cells fall under 1, they're dead)
    dP(1) = (2*a_c*s-1)*p_c*P(1)-d*P(1);
    dP(2) = 2*(1-a_c*s)*p_c*P(1)-d_c*P(2);
    if P(1)<1
        dP(1)=0;
    end
%   Leukemic clones: (if mitotic cells fall under 1, they're dead)
    for i=1:N_clones
        a_l = parameters_l(1,i);
        p_l = parameters_l(2,i);
        d_l = parameters_l(3,i);
        dP(3*i)   = (2*a_l*s-1)*p_l*P(3*i)-d*P(3*i);
        dP(3*i+1) = 2*(1-a_l*s)*p_l*P(3*i)-m*P(3*i+1)-...
                    d*P(3*i+1)-d_l*P(3*i+1);
        dP(3*i+2) = m*P(3*i+1)-d_l*P(3*i+2);
        if P(3*i)<1
            dP(3*i)=0;
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
        dP(3*i) = dP(3*i)-k_chemo(2)*p_l*P(3*i);
    end
end
