function set_publication_style(ax)
%SET_PUBLICATION_STYLE Applies a standardised, publication-quality style to an axes object.
%   This function configures font sizes, line widths, interpreters, and other
%   graphical properties to ensure a consistent and professional appearance
%   across all generated figures.
%
%   Syntax:
%       set_publication_style(ax)
%
%   Input:
%       ax - The axes handle to be styled (e.g., from gca).

    arguments
        ax (1,1) matlab.graphics.axis.Axes
    end

    set(ax, ...
        'FontSize', 22, ...                         % Set main font size for labels
        'TickLabelInterpreter', 'latex', ...          % Render tick labels using LaTeX
        'XMinorTick', 'on', ...                       % Enable minor ticks on x-axis
        'YMinorTick', 'on', ...                       % Enable minor ticks on y-axis
        'LineWidth', 1.5, ...                       % Set thickness of axis lines and ticks
        'Box', 'on', ...                            % Ensure the plot box is drawn
        'FontName', 'Times New Roman');             % Use a standard serif font

    if ~isempty(ax.Title), ax.Title.Interpreter = 'latex'; end
    if ~isempty(ax.XLabel), ax.XLabel.Interpreter = 'latex'; end
    if ~isempty(ax.YLabel), ax.YLabel.Interpreter = 'latex'; end
    
    if ~isempty(ax.Legend)
        ax.Legend.Interpreter = 'latex';
        ax.Legend.FontSize = 20;
        ax.Legend.Box = 'off';
    end
end