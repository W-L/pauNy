import plotnine as p9
import pandas as pd

# custom theme that increases font size and replaces transparent background
custom_theme = (p9.theme_minimal(base_size=16) +
                p9.theme(plot_background=p9.element_rect(fill="white", color="white"),
                         legend_title=p9.element_blank(),
                         ))


def plot_nx(nx_frame: pd.DataFrame) -> p9.ggplot:
    """
    Create an Nx curve (i.e. line plot).
    Different assemblies are mapped to color, reference boolean marker is mapped to linetype.
    Uses some standard color brewer palette for now.

    :param nx_frame: Input frame with Nx values for each assembly
    :return: plot
    """
    p = (p9.ggplot() +
         p9.geom_line(
             data=nx_frame,
             mapping=p9.aes(x='Nx', y='val', linetype="reference", color="name"),
             size=2) +
         p9.scale_color_brewer(type="qual", palette=2) +
         p9.guides(linetype=None) +
         p9.ylab("contig length (bp)") +
         p9.xlab("Nx") +
         custom_theme)
    return p




def plot_aun(aun_frame: pd.DataFrame) -> p9.ggplot:
    """
    Create a lollipop plot for the auN values.
    Draws segments and points, with assemblies mapped to color, reference boolean mapped to linetype.

    :param aun_frame: Input frame with auN values for each assembly
    :return: plot
    """
    q = (p9.ggplot(
            data=aun_frame,
            mapping=p9.aes(x='name', y='auN')) +
         p9.geom_segment(
             mapping=p9.aes(x='name', xend='name',
                            y=0, yend='auN',
                            linetype="reference",
                            color="name")) +
         p9.geom_point(
             mapping=p9.aes(color="name"),
             size=4,
             alpha=0.6) +
         p9.scale_color_brewer(type="qual", palette=2) +
         p9.guides(linetype=None, color=None) +
         p9.ylab("area under Nx") +
         p9.xlab("") +
         p9.coord_flip() +
         custom_theme)
    return q


def save_plot(plot: p9.ggplot, type: str, out_name: str, format: str = 'pdf') -> None:
    """
    Save a plot to file

    :param plot: Input plot
    :param type: Just a string that is part of the filename
    :param out_name: Base name of the output file
    :param format: Specify either pdf, png, ..
    :return: None
    """
    plot.save(f'{out_name}.{type}.{format}', dpi=600) # width=w, height=h)



