import numpy as np
import matplotlib.pylab as plt
from freeenergyframework import stats, plotlying
import itertools

def _master_plot(x, y, title='',
                 xerr=None, yerr=None,
                 method_name='', target_name='',plot_type=rf'$\Delta \Delta$ G',
                 guidelines=True, origins=True,
                 statistics=['RMSE',  'MUE'], filename=None):
    """ Handles the aesthetics of the plots in one place.

    Parameters
    ----------
    x : list
        Values to plot on the x axis
    y : list
        Values to plot on the y axis
    title : string, default = ''
        Title for the plot
    xerr : list , default = None
        Error bars for x values
    yerr : list , default = None
        Error bars for y values
    method_name : string, optional
        name of method associated with results, e.g. 'perses'
    target_name : string, optional
        name of system for results, e.g. 'Thrombin'
    plot_type : string, default = rf'$\Delta \Delta$ G'
        String of data being plotted for axis labels
    guidelines : bool, default = True
        toggles plotting of grey 0.5 and 1 kcal/mol error zone
    origins : bool, default = True
        toggles plotting of x and y axis
    statistics : list(str), default = ['RMSE',  'MUE']
        list of statistics to calculate and report on the plot
    filename : str, default = None
        filename for plot

    Returns
    -------

    """
    nsamples = len(x)
    # aesthetics
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12
    plt.rcParams['font.size'] = 12

    plt.figure(figsize=(6, 6))
    plt.subplots_adjust(left=0.2, right=0.8, bottom=0.2, top=0.8)

    plt.xlabel(f'Experimental {plot_type} ' + r'$[\mathrm{kcal\,mol^{-1}}]$')
    plt.ylabel(f'Calculated {plot_type} {method_name} ' + r'$[\mathrm{kcal\,mol^{-1}}]$')

    ax_min = min(min(x), min(y)) - 0.5
    ax_max = max(max(x), max(y)) + 0.5
    scale = [ax_min, ax_max]

    plt.xlim(scale)
    plt.ylim(scale)

    # plots x-axis and y-axis
    if origins:
        plt.plot([0, 0], scale, 'gray')
        plt.plot(scale, [0, 0], 'gray')

    #plots x=y line
    plt.plot(scale, scale, 'k:')
    if guidelines:
        small_dist = 0.5
        # plots grey region around x=y line
        plt.fill_between(scale, [ax_min - small_dist, ax_max - small_dist],
                         [ax_min + small_dist, ax_max + small_dist],
                         color='grey', alpha=0.2)
        plt.fill_between(scale, [ax_min - small_dist * 2, ax_max - small_dist * 2],
                         [ax_min + small_dist * 2, ax_max + small_dist * 2],
                         color='grey', alpha=0.2)
    # actual plotting
    cm = plt.get_cmap('coolwarm')
    clr = np.abs(x-y)
    # 2.372 kcal / mol = 4 RT
    clr = cm(clr / 2.372)
    plt.errorbar(x, y, xerr=xerr, yerr=yerr, color='gray', linewidth=0., elinewidth=2., zorder=1)
    plt.scatter(x, y,color=clr, s=10, marker='o', zorder=2)

    # stats and title
    statistics_string = ''
    for statistic in statistics:
        s = stats.bootstrap_statistic(x, y, xerr, yerr, statistic=statistic)
        string = f"{statistic}:   {s['mle']:.2f} [95%: {s['low']:.2f}, {s['high']:.2f}] " + r"$\mathrm{kcal\,mol^{-^1}}$" + "\n"
        statistics_string += string

    long_title = f'{title} \n {target_name} (N = {nsamples}) \n {statistics_string}'

    plt.title(long_title, fontsize=12, loc='right', horizontalalignment='right', family='monospace')

    if filename is None:
        plt.show()
    else:
        plt.savefig(filename,bbox_inches='tight')
    return


def plot_DDGs(results, smiles=None, method_name_x='Experimental', method_name_y='', target_name='', title='',
              map_positive=False, filename=None, symmetrise=False, plotly=False):
    """ Function to plot relative free energies

    Parameters
    ----------
    results : nx.DiGraph
        graph object with relative free energy edges
    method_name : string, optional
        name of method associated with results, e.g. 'perses'
    target_name : string, optional
        name of system for results, e.g. 'Thrombin'
    title : string, default = ''
        Title for the plot
    map_positive : bool, default=False
        whether to map all DDGs to the positive x values.
        this is an aesthetic choice
    filename : str, default = None
        filename for plot
    symmetrise : bool, default = False
        whether to plot each datapoint twice, both
        positive and negative

    Returns
    -------

    """
    # data
    if not map_positive:
        x_data = np.asarray([x.exp_DDG for x in results])
        y_data = np.asarray([x.calc_DDG for x in results])
    elif symmetrise:
        x_data = np.asarray([x.exp_DDG for x in results]+[-x.exp_DDG for x in results])
        y_data = np.asarray([x.calc_DDG for x in results]+[-x.calc_DDG for x in results])
    else:
        x_data = []
        y_data = []
        for i,j in zip([x.exp_DDG for x in results],[x.calc_DDG for x in results]):
            if i < 0:
                x_data.append(-i)
                y_data.append(-j)
            else:
                x_data.append(i)
                y_data.append(j)
        x_data = np.asarray(x_data)
        y_data = np.asarray(y_data)
    xerr = np.asarray([x.dexp_DDG for x in results])
    yerr = np.asarray([x.dcalc_DDG for x in results])

    c = np.asarray([x.other_dDDG for x in results])
    names = [f'{x.ligandA}, {x.ligandB}' for x in results]

    if plotly:
        plotlying._master_plot(x_data, y_data, c=c,
                     xerr=xerr, yerr=yerr, names=names, smiles=smiles, filename=filename, plot_type='ΔΔG',
                     title=title, method_name_x=method_name_x, method_name_y=method_name_y, target_name=target_name)
    else:
        _master_plot(x_data, y_data,
                 xerr=xerr, yerr=yerr, filename=filename,
                 title=title, method_name=method_name_y, target_name=target_name)

    return



def plot_DGs(graph, smiles=None, method_name_x='Experimental', method_name_y='', target_name='', title='', filename=None, plotly=False):
    """Function to plot absolute free energies.

    Parameters
    ----------
    graph : nx.DiGraph
        graph object with relative free energy edges
    method_name : string, optional
        name of method associated with results, e.g. 'perses'
    target_name : string, optional
        name of system for results, e.g. 'Thrombin'
    title : string, default = ''
        Title for the plot
    filename : str, default = None
        filename for plot

    Returns
    -------

    """

    # data
    x_data = np.asarray([node[1]['f_i_exp'] for node in graph.nodes(data=True)])
    y_data = np.asarray([node[1]['f_i_calc'] for node in graph.nodes(data=True)])
    xerr = np.asarray([node[1]['df_i_exp'] for node in graph.nodes(data=True)])
    yerr = np.asarray([node[1]['df_i_calc'] for node in graph.nodes(data=True)])

    # centralising
    # this should be replaced by providing one experimental result
    x_data = x_data - np.mean(x_data)
    y_data = y_data - np.mean(y_data)

    if plotly:
        plotlying._master_plot(x_data, y_data,
                 xerr=xerr, yerr=yerr,
                 origins=False, statistics=['RMSE','MUE','R2','rho'],plot_type='ΔG',
                 title=title, method_name_x=method_name_x, method_name_y=method_name_y, target_name=target_name, filename=filename)
    else:
        _master_plot(x_data, y_data,
                                   xerr=xerr, yerr=yerr,
                                   origins=False, statistics=['RMSE', 'MUE', 'R2', 'rho'], plot_type=rf'$\Delta$ G',
                                   title=title, method_name=method_name_y, target_name=target_name, filename=filename)

    return


def plot_all_DDGs(graph, smiles = None, method_name_x='Experimental', method_name_y='', target_name='', title='', filename=None, plotly=True):
    """Plots relative free energies between all ligands, which is calculated from
    the differences between all the absolute free energies. This data is different to `plot_DGs`

    Parameters
    ----------
    graph : nx.DiGraph
        graph object with relative free energy edges
    method_name : string, optional
        name of method associated with results, e.g. 'perses'
    target_name : string, optional
        name of system for results, e.g. 'Thrombin'
    title : string, default = ''
        Title for the plot
    filename : str, default = None
        filename for plot
    plotly : bool, default = True
        whether to use plotly for the plotting

    Returns
    -------

    """

    x_abs = np.asarray([node[1]['f_i_exp'] for node in graph.nodes(data=True)])
    y_abs = np.asarray([node[1]['f_i_calc'] for node in graph.nodes(data=True)])
    xabserr = np.asarray([node[1]['df_i_exp'] for node in graph.nodes(data=True)])
    yabserr = np.asarray([node[1]['df_i_calc'] for node in graph.nodes(data=True)])
    # do all to plot_all
    x_data = []
    y_data = []
    xerr = []
    yerr = []
    for a, b in itertools.combinations(range(len(x_abs)), 2):
        x = x_abs[a] - x_abs[b]
        x_data.append(x)
        x_data.append(-x)
        err = (xabserr[a]**2 + xabserr[b]**2)**0.5
        xerr.append(err)
        xerr.append(err)
        y = y_abs[a] - y_abs[b]
        y_data.append(y)
        y_data.append(-y)
        err = (yabserr[a]**2 + yabserr[b]**2)**0.5
        yerr.append(err)
        yerr.append(err)
    x_data = np.asarray(x_data)
    y_data = np.asarray(y_data)

    if plotly:
        plotlying._master_plot(x_data, y_data,
                 xerr=xerr, yerr=yerr,
                 title=title, method_name_x=method_name_x, method_name_y=method_name_y,  plot_type='ΔΔG',
                 filename=filename, target_name=target_name)

    else:
        _master_plot(x_data, y_data,
                     xerr=xerr, yerr=yerr,
                     title=title, method_name=method_name_y,
                     filename=filename, target_name=target_name)
