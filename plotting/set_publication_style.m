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

    set(ax, ...
        'FontSize', 24, ...                         % Set main font size for labels
        'TickLabelInterpreter', 'latex', ...          % Render tick labels using LaTeX
        'XMinorTick', 'on', ...                       % Enable minor ticks on x-axis
        'YMinorTick', 'on', ...                       % Enable minor ticks on y-axis
        'LineWidth', 1.5, ...                       % Set thickness of axis lines and ticks
        'Box', 'on', ...                            % Ensure the plot box is drawn
        'FontName', 'Times New Roman');             % Use a standard serif font

    % Style the title and labels
    ax.Title.Interpreter = 'latex';
    ax.XLabel.Interpreter = 'latex';
    ax.YLabel.Interpreter = 'latex';
    
    % Style the legend, if it exists
    if ~isempty(ax.Legend)
        ax.Legend.Interpreter = 'latex';
        ax.Legend.FontSize = 22;
        box(ax.Legend, 'off');
    end
end