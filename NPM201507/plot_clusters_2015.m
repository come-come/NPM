function plot_clusters_2015(D,Idx,kkk)

for c = 1:kkk
    point_c{c} = [];
end

for i = 1:length(D)
    for c = 1:kkk
        if Idx(i) == c
            point_c{c} = [point_c{c}, D{i}];
        end
    end
end

figure;

plot(point_c{1}(1,:),point_c{1}(2,:),'r.','MarkerSize',7);
hold on

if (kkk>1)&&(~isempty(point_c{2}))
    plot(point_c{2}(1,:),point_c{2}(2,:),'g.','MarkerSize',7);
    hold on
end

if kkk > 2 
    plot(point_c{3}(1,:),point_c{3}(2,:),'b.','MarkerSize',7);
    hold on
end
if kkk > 3
    plot(point_c{4}(1,:),point_c{4}(2,:),'y.','MarkerSize',7);
    hold on
end
if kkk > 4
    plot(point_c{5}(1,:),point_c{5}(2,:),'m.','MarkerSize',7);
    hold on
end
if kkk > 5
    plot(point_c{6}(1,:),point_c{6}(2,:),'k.','MarkerSize',7);
end

xlabel('Measurement 1');
ylabel('Measurement 2');

% draw a new 3D figure if there are three dimensions
if size(D{1},1) == 3
    figure;

    scatter3(point_c{1}(1,:),point_c{1}(2,:), point_c{1}(3,:), 'r.');
    hold on
    
    if (kkk>1)&&(~isempty(point_c{2}))
        scatter3(point_c{2}(1,:),point_c{2}(2,:), point_c{2}(3,:), 'g.');
        hold on
    end
    
    if kkk > 2 
        scatter3(point_c{3}(1,:),point_c{3}(2,:), point_c{3}(3,:), 'b.');
        hold on
    end
    if kkk > 3
        scatter3(point_c{4}(1,:),point_c{4}(2,:), point_c{4}(3,:), 'y.');
        hold on
    end
    if kkk > 4
        scatter3(point_c{5}(1,:),point_c{5}(2,:), point_c{5}(3,:), 'm.');
        hold on
    end
    if kkk > 5
        scatter3(point_c{6}(1,:),point_c{6}(2,:), point_c{6}(3,:), 'k.');
    end

    xlabel('Measurement 1');
    ylabel('Measurement 2');
    zlabel('Measurement 3');
end



% if ms(1) == 1
%     xlabel('Phi2');
% elseif ms(1) == 2
%     xlabel('qE');
% elseif ms(1) == 3
%     xlabel('qI');
% end
% 
% if ms(2) == 1
%     ylabel('Phi2');
% elseif ms(2) == 2
%     ylabel('qE');
% elseif ms(2) == 3
%     ylabel('qI');
% end

% if (ms(1) == 1) && (ms(2) == 2)
%     axis([-200,800,-1,9]);
% end
% if (ms(1) == 1) && (ms(2) == 3)
%     axis([-200,800,-2,8]);
% end

% %--------------
% figure;
% plot(point_c{1}(1,:),point_c{1}(2,:),'r.','MarkerSize',7);
% hold on;
% 
% if ms(1) == 1
%     xlabel('Phi2');
% elseif ms(1) == 2
%     xlabel('qE');
% elseif ms(1) == 3
%     xlabel('qI');
% end
% 
% if ms(2) == 1
%     ylabel('Phi2');
% elseif ms(2) == 2
%     ylabel('qE');
% elseif ms(2) == 3
%     ylabel('qI');
% end
% 
% % if (ms(1) == 1) && (ms(2) == 2)
% %     axis([-200,800,-1,9]);
% % end
% % if (ms(1) == 1) && (ms(2) == 3)
% %     axis([-200,800,-2,8]);
% % end
% 
% 
% figure;
% plot(point_c{2}(1,:),point_c{2}(2,:),'g.','MarkerSize',7);
% hold on;
% 
% if ms(1) == 1
%     xlabel('Phi2');
% elseif ms(1) == 2
%     xlabel('qE');
% elseif ms(1) == 3
%     xlabel('qI');
% end
% 
% if ms(2) == 1
%     ylabel('Phi2');
% elseif ms(2) == 2
%     ylabel('qE');
% elseif ms(2) == 3
%     ylabel('qI');
% end
% 
% % if (ms(1) == 1) && (ms(2) == 2)
% %     axis([-200,800,-1,9]);
% % end
% % if (ms(1) == 1) && (ms(2) == 3)
% %     axis([-200,800,-2,8]);
% % end
% 
% 
% figure;
% plot(point_c{3}(1,:),point_c{3}(2,:),'b.','MarkerSize',7);
% hold on;
% 
% if ms(1) == 1
%     xlabel('Phi2');
% elseif ms(1) == 2
%     xlabel('qE');
% elseif ms(1) == 3
%     xlabel('qI');
% end
% 
% if ms(2) == 1
%     ylabel('Phi2');
% elseif ms(2) == 2
%     ylabel('qE');
% elseif ms(2) == 3
%     ylabel('qI');
% end
% 
% % if (ms(1) == 1) && (ms(2) == 2)
% %     axis([-200,800,-1,9]);
% % end
% % if (ms(1) == 1) && (ms(2) == 3)
% %     axis([-200,800,-2,8]);
% % end
% 
% %if (ms(1) == 1) && (ms(2) == 3)
% %    axis([0,450,0,4.5]);
% %end
% 
% % if dosave == 1 
% %     fdir = [firr,'_m',num2str(ms(1)),num2str(ms(2)),'_k',num2str(kk),'_c',num2str(kkk(1)),num2str(kkk(2)),num2str(kkk(3)),'_day',num2str(cond),'.fig'];
% %     saveas(gcf,fdir);
% % end
