import numpy
import scipy
import scipy.stats
from scipy.optimize import curve_fit
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1 import make_axes_locatable


def get_cmap(N, cmap='jet'):

    """
    Returns a function that maps each index in 0, 1, ... N-1 to a distinct 
    RGB color.
    """

    import matplotlib.cm as cmx
    import matplotlib.colors as colors

    color_norm = colors.Normalize(vmin=0, vmax=N - 1)
    scalar_map = cmx.ScalarMappable(norm=color_norm, cmap=cmap)

    def map_index_to_rgb_color(index):
        return scalar_map.to_rgba(index)

    return map_index_to_rgb_color


def plot_network(G, title, filename, use_node_weights=False, use_edge_weights=False,
                 use_labels=False, use_title=False, default_size=300,
                 default_edge_thickness = 1.0, scaling_factor=1.0,
                 shell1=None, shell2=None):

    """
    
    Uses NetworkX's plotting functions to generate a figure of G.
    
    :param G: NetworkX graph object (can be directed or undirected)
    :param title: Title to apply to figure (if use_title = True)
    :param filename: Output file name
    :param use_node_weights: If true, uses a 'weight' field in G's attribute dictionary to scale node size,
    otherwise size is set to default specified in default_size parameter
    :param use_edge_weights: If true, uses a 'weight' field in G's attribute dictionary to scale edge thickness,
    otherwise thickness is set to default specified in default_edge_thickness parameter
    :param use_labels: If true, labels nodes
    :param use_title: If true, writes title into figure
    :param default_size: Default size for a node (int)
    :param default_edge_thickness: Default thickness for an edge if use_edge_weights is true
    :param scaling_factor: Multiply edge and node size parameters by X (float)
    :param shell1: If not None, then nodes specified in shell1 are the outer shell in an arranged graph (shell2 must
    not be None as well)
    :param shell2: If not None, then nodes specified in shell1 are the inner shell in an arranged graph (shell1 must
    not be None as well)
    :return: 
    """

    import networkx

    k = 1.5 / numpy.sqrt(len(G.nodes()))

    scale = 10

    if shell1 is None or shell2 is None:
        pos = networkx.spring_layout(G, k=k, scale=scale)
    else:
        pos = networkx.shell_layout(G, [shell1, shell2], scale=scale)

    node_color_test = []

    node_size_test = []

    for node in G.nodes():

        if 'color' in G.node[node]:
            # node_colors[node] = G.node[node]['color']
            node_color_test.append(G.node[node]['color'])
        else:
            # default color is red (in pajek, matplotlib, etc)
            node_color_test.append('r')

        if use_node_weights:
            node_size_test.append(abs(G.node[node].get('weight', default_size)) * scaling_factor)
        else:
            node_size_test.append(default_size)

    if use_edge_weights:

        edge_weights = []

        for edge in G.edges():
            edge_weights.append(G[edge[0]][edge[1]].get('weight', default_edge_thickness) * scaling_factor)

        networkx.draw_networkx_edges(G, pos=pos, edges=G.edges(), width=edge_weights)

    else:
        networkx.draw_networkx_edges(G, pos=pos, width=0.5)

    networkx.draw_networkx_nodes(G, pos=pos, node_color=node_color_test, node_size=node_size_test)

    if use_labels:

        labels = {}

        for node in G.nodes():
            labels[node] = node

        networkx.draw_networkx_labels(G, pos=pos, labels=labels)

    if use_title:
        plt.title(title)

    plt.xticks([])
    plt.yticks([])

    plt.tight_layout()

    plt.savefig(filename)
    plt.close()


def generate_stacked(frequency, category_breakdown,
                     xlabel, ylabel, filename, top=9, name_mapper=None, add_other=True, cmap='viridis'):

    """
    
    Generates a stacked bar plot showing relative proportions of elements described by category_breakdown sorted
    according to the frequency described in frequency dict (str : int). Helper function for stackedbar, meant to be
    called by user.
    
    :param frequency: dict (str : int) to sort keys in category_breakdown by
    :param category_breakdown: dict (str : dict : (str : int)
    :param xlabel: str, label for x-axis
    :param ylabel: str, label for y-axis
    :param filename: Output filename
    :param top: Number of entries to include (top 9 by default)
    :param name_mapper: dict (str:str) to convert categories into more user-friendly output names 
    :param add_other: If true, adds 'OTHER' category to output plot
    :param cmap: str, matplotlib to use in plot (check online for possible values)
    :return: 
    """

    freq_array = []

    if name_mapper is None:
        name_mapper = {}

    for key in frequency:
        freq_array.append((key, frequency[key]))

    # sort to get descending values
    freq_array = sorted(freq_array, key=lambda freq_array: freq_array[1], reverse=True)

    # top X entries
    if top != -1:
        top_types = [d[0] for d in freq_array[0:top - 1]]
    else:
        top_types = [d[0] for d in freq_array]

    # make a stacked bar chart with products

    keys = category_breakdown.keys()
    bar_data = []
    key_order = sorted(keys)
    type_set = set()

    for key in key_order:

        data = defaultdict(int)

        for ptype in category_breakdown[key]:
            if ptype in top_types:
                data[ptype] += 1
                type_set.add(ptype)
            else:
                data['OTHER'] += 1

        value_array = []

        for ptype in data:
            value_array.append(data[ptype])

        sum_types = numpy.sum(value_array)

        value_array = []

        for ptype in data:
            data[ptype] = float(data[ptype]) / float(sum_types)
            value_array.append(data[ptype])

        bar_data.append(data)

    order = sorted([item for item in list(type_set)])

    if add_other and 'OTHER' not in order:
        order.append('OTHER')

    stackedbar(order,
               bar_data,
               key_order,
               xlabel,
               ylabel,
               filename,
               rotation='vertical',
               mapper=name_mapper,
               cmap=cmap)


def stackedbar(value_keys, dict_array, xtick_labels, xlabel, ylabel, filename, rotation='horizontal', mapper=None,
               cmap='jet'):

    """
    
    Generates a stacked bar plot using the output of generate_stacked (usually).
    
    :param value_keys: ordered list to get values from each dict_array
    :param dict_array: ordered list of dicts
    :param xtick_labels: labels for x-axis
    :param xlabel: x-axis label
    :param ylabel: y-axis label
    :param filename: output file name
    :param rotation: rotation of x-axis labels
    :param mapper: Convert xtick_labels into more user-friendly labels
    :param cmap: Matplotlib cmap to use (search online for options)
    :return: 
    """

    color_mapper = get_cmap(len(value_keys), cmap=cmap)

    if mapper is None:
        mapper = {}

    try:

        plt.subplots()

        index = numpy.arange(len(xtick_labels))
        bar_width = 0.55

        opacity = 0.8
        previous_top = []
        counter = 0

        patch_array = []

        for key in value_keys:

            data_array = []

            for data_dict in dict_array:
                # extract the data in the specified order
                data_array.append(data_dict[key])

            if len(previous_top) == 0:
                previous_top = [0] * len(data_array)

            plt.bar(index, data_array, bar_width, alpha=opacity, color=color_mapper(counter), align='center',
                    bottom=previous_top)

            # need to keep track of where we are on the plot.
            for i in range(0, len(data_array)):
                previous_top[i] = previous_top[i] + data_array[i]

            artist = mpatches.Patch(color=color_mapper(counter), label=mapper.get(key, str(key).upper()))
            patch_array.append(artist)

            counter = counter + 1

        plt.axis('tight')

        # this usually looks okay but can be tweaked to be better
        plt.legend(handles=patch_array, loc='lower center',
                   bbox_to_anchor=(0.5, min(-0.036 * float(len(value_keys)), -0.3)), ncol=3)

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        plt.xticks(index, [mapper.get(k, k) for k in xtick_labels], rotation=rotation)

        plt.savefig(filename, bbox_inches='tight')
        plt.close('all')
    except:
        plt.close('all')
        raise


def generate_heatmap(matrix, names, ylabel, filename, vmin=None, vmax=None, mapper=None, use_labels=True, cmap='jet',
                     plot_color_bar=True, x_label_rotation='vertical', y_label_rotation='horizontal'):

    """
    
    Generates a symmetric heatmap using the data in matrix (NxN matrix) using labels obtained from names.
    
    :param matrix: NxN matrix of floats
    :param names: Nx1 column and row names to display
    :param ylabel: Label to the right of the added colorbar
    :param filename: Output file name
    :param vmin: Minimum value for the colorbar
    :param vmax: Maximum value for the colorbar
    :param mapper: Converts names to user-friendly values if present in mapper
    :param use_labels: If true, labels columns and rows in the plot
    :param cmap: Matplotlib colormap
    :param plot_color_bar: If true, plots colorbar.
    :param x_label_rotation: str, 'horizontal' or 'vertical'
    :param y_label_rotation: str, 'horizontal' or 'vertical'
    :return: 
    """

    generate_heatmap_asymmetric(matrix,
                                names,
                                names,
                                ylabel,
                                filename,
                                vmin,
                                vmax,
                                mapper,
                                cmap=cmap,
                                plot_color_bar=plot_color_bar,
                                dim_override=None,
                                patches=None,
                                x_label_rotation=x_label_rotation,
                                y_label_rotation=y_label_rotation)


def generate_heatmap_asymmetric(matrix, x_names, y_names, ylabel, filename, vmin=None, vmax=None, mapper=None,
                                cmap='viridis', plot_color_bar=True, dim_override=None,
                                patches=None, x_label_rotation = 'horizontal', y_label_rotation='horizontal'):

    """
    
    Similar to generate_heatmap, but does not assume matrix is symmetrical.
    
    :param matrix: NxM matrix of floats
    :param x_names: Mx1 list to label columns
    :param y_names: Nx1 list to label rows
    :param ylabel: Label to the right of the added colorbar
    :param filename: Output file name
    :param vmin: Minimum value for the colorbar
    :param vmax: Maximum value for the colorbar
    :param mapper: Converts names to user-friendly values if present in mapper
    :param cmap: Matplotlib colormap
    :param plot_color_bar: If true, plots colorbar.
    :param dim_override: If not none, assumes dim_override is a tuple (float, float) and fixes the canvas size to match
    that input.
    :param patches: If not none, assumes patches is a list of (row, col, _, color) tuples and draws a star on the 
    heatmap at those coordinates
    :param x_label_rotation: str, 'horizontal' or 'vertical'
    :param y_label_rotation: str, 'horizontal' or 'vertical'
    :return: 
    """

    try:
        fig, ax = plt.subplots()
        if vmin is not None and vmax is not None:

            if patches is None:
                plt.pcolormesh(matrix, zorder=1, vmin=vmin, vmax=vmax, cmap=cmap,  edgecolors='k')
            else:
                plt.pcolormesh(matrix, zorder=1, vmin=vmin, vmax=vmax, cmap=cmap, edgecolors='k')
                for (r, c, _, color) in patches:
                    plt.text(c + 0.5,
                             r + 0.25,
                             '*',
                             horizontalalignment='center',
                             verticalalignment='center',
                             color=color,
                             fontweight='bold')
        else:
            if patches is None:
                plt.pcolormesh(matrix, zorder=1, cmap=cmap, edgecolors='k')
            else:
                plt.pcolormesh(matrix, zorder=1, cmap=cmap, edgecolors='k')
                for (r, c, _, color) in patches:
                    plt.text(c + 0.5,
                             r + 0.25,
                             '*',
                             horizontalalignment='center',
                             verticalalignment='center',
                             color=color,
                             fontweight='bold')

        if dim_override is not None:
            fig = plt.gcf()
            fig.set_size_inches(w=dim_override[0], h=dim_override[1])

        plt.axis('tight')
        if mapper is not None:
            x_names = [mapper.get(name, name) for name in x_names]
            y_names = [mapper.get(name, name) for name in y_names]
        else:
            x_names = list(x_names)
            y_names = list(y_names)

        plt.xticks(numpy.arange(len(x_names)) + 0.5, x_names, rotation='horizontal')
        plt.yticks(numpy.arange(len(y_names)) + 0.5, y_names, rotation='vertical')
        ax.set_yticklabels(y_names, rotation=y_label_rotation, va='center')
        ax.set_xticklabels(x_names, rotation=x_label_rotation)

        if plot_color_bar:

            ax = plt.gca()
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('right', size=0.2, pad=0.1)
            cb = plt.colorbar(cax=cax)
            cb.set_label(ylabel)

        plt.savefig(filename, bbox_inches='tight')
        plt.close('all')

    except:
        plt.close('all')
        raise


def lineplot(xvalues, yvalues, legends, xlabel, ylabel, filename, yerrors=None, ymin=0,
             showx=True):

    """
    
    Generates a line plot
    
    :param xvalues: x_values to plot
    :param yvalues: y_values to plot
    :param legends: list of strings to label line plots
    :param xlabel: x label
    :param ylabel: y label
    :param filename: Output file name
    :param yerrors: If not none, a list equal in length to yvalues to plot error bars
    :param ymin: minimum value for y on the plot
    :param showx: If true, shows x tick labels
    :return: 
    """

    # yvalues is a list of lists, legends = len(yvalues), xvalues = len(yvalues[i])
    # if provided, yerrors also == len(yvalues)

    colors = get_cmap(len(yvalues))

    try:

        fig, ax = plt.subplots()

        y_max = 0

        frame1 = plt.gca()

        if yerrors is None:
            for y, i in zip(yvalues, range(0, len(yvalues))):
                ax.plot(xvalues, y, color=colors(i), label=legends[i])
                if max(y) > y_max:
                    y_max = max(y)
        else:
            for y, i in zip(yvalues, range(0, len(yvalues))):
                plt.errorbar(xvalues, y, yerr=yerrors[i], color=colors(i), label=legends[i])
                if max(y) > y_max:
                    y_max = max(y)

        ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)

        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc='best')

        frame1.axes.get_xaxis().set_visible(showx)

        print [ymin, y_max]

        ax.set_ylim([ymin, y_max * 1.1])

        plt.savefig(filename, bbox_inches='tight')
        plt.close('all')

    except:

        plt.close('all')
        raise


def linebar(keys, values, yvalues, xlabel, y1label, y2label, filename, mapper=None, edgecolor='black',
            show_ticks=True):

    """
    
    Plots an overlapping bar and line plot.
    
    :param keys: x axis keys
    :param values: values for bar chart
    :param yvalues: values for line plot
    :param xlabel: Label for x-axis
    :param y1label: Label for left y-axis (bar chart)
    :param y2label: Label for right y-axis (line plot)
    :param filename: Output file name
    :param mapper: Dict (str: str) to convert x-axis labels to more user-friendly ones
    :param edgecolor: Edge color for bars (black, red, etc...)
    :param show_ticks: If true, shows labels on x/y-axis plots
    :return: 
    """

    try:

        fig, ax = plt.subplots()

        index = numpy.arange(len(keys))
        bar_width = 0.7

        opacity = 0.8

        ax.bar(index, values, bar_width,
               alpha=opacity,
               color='b',
               align='center',
               edgecolor=edgecolor)

        ax.set_ylabel(y1label)
        ax.set_xlabel(xlabel)

        ax.set_xticks(index)

        if not show_ticks:
            plt.tick_params(
                axis='x',  # changes apply to the x-axis
                which='both',  # both major and minor ticks are affected
                bottom='off',  # ticks along the bottom edge are off
                top='off',  # ticks along the top edge are off
                labelbottom='off')  # labels along the bottom edge are off

        ax2 = ax.twinx()
        ax2.plot(index, yvalues, 'r-', linewidth=2.0)
        ax2.set_ylabel(y2label, color='red')

        # turns off offset on graph if change between points is small
        ax2.get_yaxis().get_major_formatter().set_useOffset(False)

        ax2.yaxis.label.set_color('red')
        ax2.tick_params(axis='y', colors='red')

        if mapper is None:
            ax.set_xticklabels([str(key) for key in keys], rotation='vertical')
        else:
            ax.set_xticklabels([str(mapper[k]) for k in keys], rotation='vertical')

        ax.set_ylim(0, max(values) * 1.1)
        ax2.set_ylim(min(yvalues) * 0.9, max(yvalues) * 1.1)

        # plt.axis('tight')
        plt.savefig(filename, bbox_inches='tight')
        plt.close('all')

    except:
        plt.close('all')
        raise


def multi_bargraph(key_list, value_list, legend, xlabel, ylabel, filename, rotation='horizontal', mapper=None):

    """
    
    Plots multiple data series as a bar chart. Bars are scaled by the length of value_list.
    
    :param key_list: x-axis labels
    :param value_list: list of list (equal of length to key_list)
    :param legend: Legend entries (can be empty list)
    :param xlabel: Label for x-axis
    :param ylabel: Label for y-axis
    :param filename: Output filename
    :param rotation: Rotation of x-axis labels
    :param mapper: dict (str: str) to convert x-axis labels
    :return: 
    """

    color_order = ['c', 'r', 'y', 'g', 'b']

    try:

        plt.subplots()

        counter = 0

        lhandles = []

        for keys, values, color, lentry in zip(key_list, value_list, color_order, legend):
            index = numpy.arange(len(keys))
            bar_width = 0.35

            opacity = 0.8
            rects1 = plt.bar(index - float(bar_width) / 2.0 + bar_width * counter, values, bar_width,
                             alpha=opacity,
                             color=color,
                             align='center',
                             label=lentry)

            lhandles.append(rects1)

            counter += 1

        plt.legend(handles=lhandles, loc='best', ncol=1)

        temp_ys = []

        for vlist in value_list:
            temp_ys.extend(vlist)

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        if mapper is None:
            plt.xticks(index, keys, rotation=rotation)
        else:
            plt.xticks(index, [mapper.get(k, k) for k in keys], rotation=rotation)

        plt.savefig(filename, bbox_inches='tight')
        plt.close('all')

        fhandle = open(filename + "_key labels.txt", 'w')
        for key in keys:
            fhandle.write(key + "\n")
        fhandle.close()

    except:
        plt.close('all')
        raise


# here values is a iterable of lists
def multi_histogram(values, xlabel, ylabel, filename, use_median_filter=False,
                    legends=None, apply_norm=False):

    """
    
    Plots a set of histograms.
    
    :param values: list of lists of data
    :param xlabel: Label of x-axis
    :param ylabel: Label for y-axis (usually 'count')
    :param filename: Output filename
    :param use_median_filter: If true, sets boundary of last bin to 5x median value of the data, otherwise bin
    boundaries are min/max of data.
    :param legends: Legend entry for each dataset (can be empty)
    :param apply_norm: If true, normalizes data such that sum of AUC is 1.
    :return: 
    """

    colors = ['b', 'r', 'g', 'y']

    try:

        plt.subplots()

        first_bins = []

        patches = []

        for i in range(0, len(values)):

            if not use_median_filter:
                vrange = [min(values[i]), max(values[i])]
            else:
                vrange = [min(values[i]), 5.0 * numpy.median(values[i])]

            if i == 0:
                hist, first_bins = numpy.histogram(values[i], range=vrange, normed=apply_norm)

            if legends is not None:
                patches.append(mpatches.Patch(color=colors[i], label=legends[i]))

            plt.hist(values[i], color=colors[i], bins=first_bins, normed=apply_norm, alpha=0.8)

        plt.axis('tight')

        plt.legend(handles=patches, loc='upper right')

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        plt.savefig(filename, bbox_inches='tight')
        plt.close('all')

    except:
        plt.close('all')
        raise


def boxplot(xvalues, yvectors, xlabel, ylabel, filename, showfliers=False):

    """
    
    Interface for matplotlib boxplot function.
    
    :param xvalues: xvalues
    :param yvectors: yvalues at each x-value
    :param xlabel: Label for x axis
    :param ylabel: Label for y axis
    :param filename: Output filename
    :param showfliers: Show outliers ('fliers')
    :return: 
    """

    try:

        plt.subplots()

        plt.boxplot(yvectors, labels=xvalues, showfliers=showfliers)

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        plt.savefig(filename, bbox_inches='tight')
        plt.close('all')

    except:
        plt.close('all')
        raise


def scatter(xvalues, yvalues, xlabel, ylabel, filename, regression=None):

    """
    
    Scatter plot with optional regressions.
    
    :param xvalues: x values
    :param yvalues: y values
    :param xlabel: Label for x axis
    :param ylabel: Label for y axis
    :param filename: Output filename
    :param regression: Can be 'linear' 'exponential' 'logistic'; if None no regression is plotted. Logistic regressions
    are using the logistic function, linear is ax+b, exponential is A*exp(b*t).
    :return: 
    """

    try:

        plt.subplots()

        plt.scatter(xvalues, yvalues, color='b')

        if regression is not None:

            lsp = numpy.linspace(min(xvalues), max(xvalues))

            if regression == 'linear':

                slope, intercept, Rsq, pvalue, stderr = scipy.stats.linregress(xvalues, yvalues)

                reg_values = [slope * pos + intercept for pos in lsp]

                equation_str = '%0.3f*t + %0.2f' % (slope, intercept)

                equation_label = 'Linear'

            elif regression == 'exponential':

                def test_func(x, a, b, c):
                    return a * numpy.exp(b * x) + c

                popt, pcov = curve_fit(test_func, numpy.asarray(xvalues), numpy.asarray(yvalues), p0=[6, 0.3, 5])

                a = popt[0]
                b = popt[1]
                c = popt[2]

                if c > 0:
                    equation_str = '%i exp(%0.3f*t) + %i' % (int(a), b, int(c))
                else:
                    equation_str = '%i exp(%0.3f*t) - %i' % (int(a), b, int(c) * -1)
                equation_label = 'Exponential'

                ybar = numpy.mean(yvalues)
                sstot = 0
                ssres = 0

                for x, y in zip(xvalues, yvalues):
                    sstot = sstot + (test_func(x, a, b, c) - ybar) ** 2
                    ssres = ssres + (test_func(x, a, b, c) - y) ** 2

                Rsq = 1 - ssres / sstot

                reg_values = [popt[0] * numpy.exp(popt[1] * x) + popt[2] for x in lsp]

            elif regression == 'logistic':

                def test_func(x, mu, s):
                    return 1 / (1 + numpy.exp(-(x - mu) / s))

                popt, pcov = curve_fit(test_func, numpy.asarray(xvalues), numpy.asarray(yvalues), p0=[5, 1])

                mu = popt[0]
                s = popt[1]

                equation_str = '1/(1+e^(-(x-%0.2f)/%0.2f)' % (mu, s)
                equation_label = 'Logistic, midpoint: %0.2f' % mu

                ybar = numpy.mean(yvalues)
                sstot = 0
                ssres = 0

                for x, y in zip(xvalues, yvalues):
                    sstot = sstot + (test_func(x, mu, s) - ybar) ** 2
                    ssres = ssres + (test_func(x, mu, s) - y) ** 2

                Rsq = 1 - ssres / sstot

                reg_values = [test_func(x, mu, s) for x in lsp]
            else:
                raise ValueError('Unknown regression type')

            plt.annotate(equation_str, (0.05, 0.9), xycoords='axes fraction')
            plt.annotate('Rsq = {0:.3f}'.format(Rsq), (0.05, 0.85), xycoords='axes fraction')
            plt.plot(lsp, reg_values, color='r', label=equation_label, linewidth=3)

        plt.axis('tight')

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        plt.savefig(filename, bbox_inches='tight')
        plt.close('all')

    except:
        plt.close('all')
        raise


def plotcdf(values, xlabel, ylabel, legend, filename, plot_survival=True):

    """
    
    Plots a cumulative distribution function for values.
    
    :param values: data
    :param xlabel: Label for x-axis
    :param ylabel: Label for y-axis
    :param legend: Legend for plot (can be empty list)
    :param filename: Output filename
    :param plot_survival: If true, plots 1-y-the survival function
    :return: 
    """

    try:

        fig, ax = plt.subplots()

        d, base = numpy.histogram(values, bins=50)
        cumulative = numpy.cumsum(d)

        plt.plot(base[:-1], cumulative, label=legend, linewidth=2.0)

        if plot_survival:
            plt.plot(base[:-1], len(values) - cumulative, label='Survival (1-y)', linewidth=2.0)

        plt.axis('tight')

        handles, labels = ax.get_legend_handles_labels()

        ax.legend(handles, [legend, 'Survival (1-y)'], loc='best')

        plt.xlim([0, 100.1])

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        plt.savefig(filename, bbox_inches='tight')
        plt.close('all')

    except:
        plt.close('all')
        raise


def histogram(values, xlabel, ylabel, filename, use_median_filter=False):

    """

    Plots a histogram.

    :param values: list of lists of data
    :param xlabel: Label of x-axis
    :param ylabel: Label for y-axis (usually 'count')
    :param filename: Output filename
    :param use_median_filter: If true, sets boundary of last bin to 5x median value of the data, otherwise bin
    boundaries are min/max of data.
    :return: 
    """

    try:

        plt.subplots()
        values_array = list(values)

        if use_median_filter:
            plt.hist(values_array, color='b', range=[min(values_array), 5.0 * numpy.median(values_array)])
        else:
            plt.hist(values_array, color='b')

        plt.axis('tight')

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        plt.savefig(filename, bbox_inches='tight')
        plt.close('all')

    except:
        plt.close('all')
        raise


def bargraph(keys, values, xlabel, ylabel, filename, rotation='horizontal', mapper=None, yerr=None):

    """
    
    Plots a standard bar chart.
    
    :param keys: labels for x-axis data points
    :param values: value at each x-axis point
    :param xlabel: Label for x-axis
    :param ylabel: Label for y-axis
    :param filename: Output filename
    :param rotation: Rotation of x-axis labels
    :param mapper: Dict (str:str) to convert labels
    :param yerr: If not None, errors to plot on each bar (length = len(keys))
    :return: 
    """

    try:

        plt.subplots()

        index = numpy.arange(len(keys))
        bar_width = 0.55

        opacity = 0.8

        plt.bar(index, values, bar_width,
                         alpha=opacity,
                         color='b',
                         align='center')

        if yerr is not None:
            plt.errorbar(index, values, yerr=yerr, fmt='none', color='k')

        plt.axis('tight')

        if yerr is not None:
            max_index = numpy.where(values == max(values))
            plt.ylim([0, max(values) + yerr[max_index[0]] * 1.2])

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        if mapper is None:
            plt.xticks(index, keys, rotation=rotation)
        else:
            plt.xticks(index, [mapper.get(k, mapper.get(k.upper(), k)) for k in keys], rotation=rotation)

        plt.savefig(filename, bbox_inches='tight')
        plt.close('all')

    except:
        plt.close('all')
        raise
