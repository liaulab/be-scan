try: 
    from .scatterplot import plot_scatterplot
    from .correlation_scatter import plot_corr_scatterplot
    from .correlation_heatmap import plot_corr_heatmap
    from .boxes import plot_boxes
except ImportError: 
    print('relative import failed')

try: 
    from scatterplot import plot_scatterplot
    from correlation_scatter import plot_corr_scatterplot
    from correlation_heatmap import plot_corr_heatmap
    from boxes import plot_boxes
except ImportError: 
    print('absolute import failed')
