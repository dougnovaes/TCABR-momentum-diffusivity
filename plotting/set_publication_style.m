function set_publication_style(ax)
%SET_PUBLICATION_STYLE Applies a standardised, publication-quality style to an axes object.
%   Configures font sizes, line widths, interpreters, and other graphical
%   properties for a consistent and professional appearance.

    arguments
        ax (1,1) matlab.graphics.axis.Axes
    end

    set(ax, ...
        'FontSize', 22, ...
        'TickLabelInterpreter', 'latex', ...
        'XMinorTick', 'on', ...
        'YMinorTick', 'on', ...
        'LineWidth', 1.5, ...
        'Box', 'on', ...
        'FontName', 'Times New Roman');

    if ~isempty(ax.Title), ax.Title.Interpreter = 'latex'; end
    if ~isempty(ax.XLabel), ax.XLabel.Interpreter = 'latex'; end
    if ~isempty(ax.YLabel), ax.YLabel.Interpreter = 'latex'; end
    
    if ~isempty(ax.Legend)
        ax.Legend.Interpreter = 'latex';
        ax.Legend.FontSize = 20;
        ax.Legend.Box = 'off';
    end
end