function x = plot_heatmap(input_data, title_, filename, min, max, xrange, yrange)

    imagesc(mean(input_data,3),[min,max])
    colorbar
    title(title_)
    xlabel("% Flipped in J, one to zero")
    ylabel("% Flipped in I, zero to one")
    xticklabels = xrange;
    xticks = linspace(1, size(mean(input_data,3), 2), numel(xticklabels));
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
    yticklabels = yrange;
    yticks = linspace(1, size(mean(input_data,3), 1), numel(yticklabels));
    set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);
    set(gca, 'XAxisLocation', 'top');
    m = round(mean(input_data,3),2);
    for i = 1:length(yrange)
        for j = 1:length(xrange)
            text(xticks(j),yticks(i),num2str(m(i,j)));
        end
    end
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(filename,'-bestfit','-dpdf')
    
    x = true;

end


%
% Original Code for this (long-form)

% imagesc(mean(es_lr,3),[0,.55])
% colorbar
% title('MEAN ES\_LR')
% xlabel("% Flipped in J, one to zero")
% ylabel("% Flipped in I, zero to one")
% xticklabels = EXP_RANGE_J;
% xticks = linspace(1, size(mean(es_lr,3), 2), numel(xticklabels));
% set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
% yticklabels = EXP_RANGE;
% yticks = linspace(1, size(mean(es_lr,3), 1), numel(yticklabels));
% set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);
% set(gca, 'XAxisLocation', 'top');
% m = round(mean(es_lr,3),2);
% for i = 1:length(EXP_RANGE)
%     for j = 1:length(EXP_RANGE_J)
%         text(xticks(j),yticks(i),num2str(m(i,j)));
%     end
% end
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print('figure1','-bestfit','-dpdf')