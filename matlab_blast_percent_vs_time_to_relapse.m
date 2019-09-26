function matlab_blast_percent_vs_time_to_relapse
    BM_blast_Ding               = [60   12   20   54   81];
    Time_to_relapse_Ding        = [252  805  505  365  235];
    BM_blast_past_Shlush        = [70   35   43   95];
    Time_to_relapse_past_Shlush = [480  180  480  240];
    Time_to_relapse_next_Shlush = [210  300];
    ID_next_Shlush              = [1    15];

    figure(1);clf

    p(1) = scatter(Time_to_relapse_Ding,BM_blast_Ding,'filled','DisplayName','Ding et al.');hold on
        p(1).SizeData           = 200;
        p(1).MarkerFaceColor    = 'r';
        p(1).LineWidth          = 2;
        p(1).MarkerEdgeColor    = 'k';
    p(2) = scatter(Time_to_relapse_past_Shlush,BM_blast_past_Shlush,'filled','DisplayName','Shlush et al. (actual data)');
        p(2).SizeData           = 200;
        p(2).MarkerFaceColor    = 'b';
        p(2).LineWidth          = 2;
        p(2).MarkerEdgeColor    = 'k';

    x_regression    = [Time_to_relapse_Ding     Time_to_relapse_past_Shlush]';
    y_regression    = [BM_blast_Ding            BM_blast_past_Shlush]';

    [Blast,Time_to_relapse,PB_blast_next_Shlush] = linear_regression(x_regression,y_regression,Time_to_relapse_next_Shlush);
    p(3) = plot(Time_to_relapse,Blast,'k--','LineWidth',2,'DisplayName','Linear regression');

    p(4) = scatter(Time_to_relapse_next_Shlush,PB_blast_next_Shlush,'filled','DisplayName','Shlush et al. (estimation)');
        p(4).SizeData           = 200;
        p(4).MarkerFaceColor    = 'g';
        p(4).LineWidth          = 2;
        p(4).MarkerEdgeColor    = 'k';

    for i=1:length(ID_next_Shlush)
        fprintf('BM blast for patient %d is %f\n',ID_next_Shlush(i),PB_blast_next_Shlush(i));
    end

    set(gca,'fontsize',20);
    ylim([0 100]);
    l = xlabel('Time to relapse (days)');
    l.FontSize  = 25;
    l   = ylabel('%BM Blast');
    l.FontSize  = 25;

    [l,objh]    = legend(p(1:4),'Location','northeast');
    l.FontSize  = 20;
    ch = findobj(objh,'type','patch');
    set(ch,'Markersize',10);

    figure_name = ['plot_blast_percent_vs_time_to_relapse_all.png'];
    saveas(gcf,figure_name);
end
function [Blast,T_relapse,y_next] = linear_regression(x_regression,y_regression,x_next)
    fun         = @(parameter_set,Time_to_relapse) linear_func(parameter_set,Time_to_relapse);
    [para,resnorm]= lsqcurvefit(fun,[10,10],x_regression,y_regression);
    T_relapse   = [0:100:900];
    A           = para(1);
    B           = para(2);
    Blast       = A + B*T_relapse;
    y_next      = A + B*x_next;
end
function Blast = linear_func(parameter_set,Time_to_relapse)
    A   = parameter_set(1);
    B   = parameter_set(2);
    Blast   = A + B*Time_to_relapse;
end
