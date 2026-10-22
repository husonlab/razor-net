#!/usr/bin/env python3
"""Draws the figures of the two simulation studies and writes their summary tables.

    figures.py <simulation-study-1.tsv> <simulation-study-2.tsv> <figure directory> <table directory>

The inputs are the tables written by collect.py. Writes one figure per study and number of taxa or sequences,
the two summary tables, and simulation-numbers.tex with the numbers that the text quotes. Needs matplotlib.
"""
import collections
import csv
import math
import os
import sys

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.ticker import FuncFormatter

INK = '#0b0b0b'
SECONDARY = '#52514e'
GRID = '#e4e3df'
REFERENCE = '#9b9a94'
# one colour and one marker per method, the same in every figure
STYLE = {
    'RazorNet': ('#2a78d6', 'o'), 'RazorHaploNet': ('#2a78d6', 'o'),
    'NeighborNet': ('#eb6834', 's'), 'MedianJoining': ('#eb6834', 's'),
    'RazorNet+Simplify': ('#1baf7a', '^'), 'RazorHaploNet+Simplify': ('#1baf7a', '^'),
    'TCS': ('#eda100', 'D'),
}
NAMES = {
    'RazorNet': 'RazorNet', 'RazorNet+Simplify': 'RazorNet + Simplify', 'NeighborNet': 'Neighbor-Net',
    'Distances': 'the distances themselves', 'RazorHaploNet': 'RazorHaploNet',
    'RazorHaploNet+Simplify': 'RazorHaploNet + Simplify', 'MedianJoining': 'Median joining', 'TCS': 'TCS',
    'Truth': 'True tree',
}

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'], 'font.size': 8,
    'axes.edgecolor': SECONDARY, 'axes.labelcolor': INK, 'axes.linewidth': 0.6, 'axes.titlesize': 8.5,
    'axes.titleweight': 'bold', 'text.color': INK, 'xtick.color': SECONDARY, 'ytick.color': SECONDARY,
    'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5, 'legend.fontsize': 8, 'pdf.fonttype': 42,
})


def number(text):
    return float(text) if text not in ('', None) else math.nan


def mean(values):
    values = [v for v in values if v == v]
    return sum(values) / len(values) if values else math.nan


def standard_error(values):
    values = [v for v in values if v == v]
    if len(values) < 2:
        return 0.0
    m = sum(values) / len(values)
    return math.sqrt(sum((v - m) ** 2 for v in values) / (len(values) - 1) / len(values))


def quantile(values, q):
    values = sorted(v for v in values if v == v)
    if not values:
        return math.nan
    position = q * (len(values) - 1)
    low = int(math.floor(position))
    high = min(low + 1, len(values) - 1)
    return values[low] + (values[high] - values[low]) * (position - low)


def median(values):
    return quantile(values, 0.5)


def read(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def ratio(row, numerator, denominator):
    d = number(row[denominator])
    return number(row[numerator]) / d if d == d and d > 0 else math.nan


def style_axes(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(length=2.5, width=0.6)


def draw_series(ax, xs, values, method, positions=None, use_median=False):
    """one line per method: the mean over the replicates with its standard error, or the median with the quartiles"""
    if use_median:
        means = [median(values[x]) for x in xs]
        errors = [[m - quantile(values[x], 0.25) for m, x in zip(means, xs)],
                  [quantile(values[x], 0.75) - m for m, x in zip(means, xs)]]
    else:
        means = [mean(values[x]) for x in xs]
        errors = [standard_error(values[x]) for x in xs]
    where = positions if positions is not None else xs
    if method in STYLE:
        colour, marker = STYLE[method]
        ax.errorbar(where, means, yerr=errors, fmt='none', ecolor=colour, elinewidth=1.0, capsize=2.5, capthick=1.0, zorder=3)
        ax.plot(where, means, color=colour, linewidth=1.5, marker=marker, markersize=6, markeredgecolor='white',
                markeredgewidth=1.2, label=NAMES[method], zorder=4)
    else:  # a reference, not a method: thin and grey
        ax.plot(where, means, color=REFERENCE, linewidth=1.0, label=NAMES[method], zorder=2)


def percent(value, _):
    return f'{100 * value:g}%'


# ---------------------------------------------------------------------------------------------------------
# study 1: conflicting gene trees
# ---------------------------------------------------------------------------------------------------------
SETTINGS_1 = [('inf', 'infinite sites'), ('low', 'low divergence'), ('mid', 'medium divergence'), ('high', 'high divergence')]
METHODS_1 = ['Distances', 'NeighborNet', 'RazorNet', 'RazorNet+Simplify']


def study1(rows, figure_directory, table_directory):
    for row in rows:
        row['n'], row['k'] = (int(x) for x in row['cell'][1:].split('k'))
        row['conflictRecall'] = ratio(row, 'sharedConflict', 'trueConflict')
        row['conflictPrecision'] = ratio(row, 'sharedConflict', 'shownConflict')
        row['recall'] = ratio(row, 'sharedSupport', 'trueSupport')
        row['precision'] = ratio(row, 'sharedSupport', 'shownSupport')
        # support displayed for splits that are not true, and support displayed for true conflicting splits beyond
        # their true support, both relative to the total true support
        row['falseSplit'] = ratio(row, 'falseSupport', 'trueSupport')
        true = number(row['trueSupport'])
        row['overDisplay'] = ((number(row['shownConflict']) - number(row['sharedConflict']) - number(row['falseSupport'])) / true
                              if true == true and true > 0 else math.nan)
    done = [r for r in rows if r['status'] == 'ok']
    sizes = sorted(set(r['n'] for r in rows))

    for n in sizes:
        figure, axes = plt.subplots(2, 4, figsize=(7.2, 3.9), sharex=True, sharey='row')
        for column, (setting, title) in enumerate(SETTINGS_1):
            for method in METHODS_1:
                recall = collections.defaultdict(list)
                false = collections.defaultdict(list)
                for r in done:
                    if r['setting'] == setting and r['n'] == n and r['method'] == method:
                        if r['k'] >= 2:
                            recall[r['k']].append(r['conflictRecall'])
                        false[r['k']].append(r['falseSplit'])
                if recall:
                    draw_series(axes[0][column], sorted(recall), recall, method)
                if false:
                    draw_series(axes[1][column], sorted(false), false, method)
            axes[0][column].set_title(title)
            axes[1][column].set_xlabel('gene trees')
            axes[1][column].set_xticks([1, 2, 3, 4, 5, 6])
            for ax in (axes[0][column], axes[1][column]):
                style_axes(ax)
        axes[0][0].set_ylabel('conflict found\n(share of the true conflict)')
        axes[0][0].set_ylim(0, 1.05)
        axes[0][0].yaxis.set_major_formatter(FuncFormatter(percent))
        axes[1][0].set_ylabel('false conflict\n(relative to the true support)')
        # logarithmic above 0.1%, linear below, so that a false conflict of zero has a place
        axes[1][0].set_yscale('symlog', linthresh=0.001, linscale=0.5)
        axes[1][0].set_ylim(-0.0003, 5)
        axes[1][0].set_yticks([0, 0.001, 0.01, 0.1, 1])
        axes[1][0].yaxis.set_major_formatter(FuncFormatter(percent))
        axes[1][0].yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        handles, labels = axes[1][0].get_legend_handles_labels()
        figure.legend(handles, labels, loc='upper center', ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.0),
                      handlelength=2.2, columnspacing=1.8)
        figure.tight_layout(rect=(0, 0, 1, 0.94))
        figure.savefig(os.path.join(figure_directory, f'sim-trees-n{n}.pdf'))
        plt.close(figure)

    # the summary table: one row per setting, number of taxa and method
    with open(os.path.join(table_directory, 'simulation-study-1-summary.tsv'), 'w') as f:
        f.write('setting\ttaxa\tmethod\tdatasets\ttimeouts\tseconds\tnodes\tcycles\tlength\tfit\trecall\tprecision'
                '\tconflictRecall\tconflictPrecision\tfalseSplit\toverDisplay\tfalseSplitOneTree\n')
        for setting, _ in SETTINGS_1:
            for n in sizes:
                for method in METHODS_1:
                    every = [r for r in rows if r['setting'] == setting and r['n'] == n and r['method'] == method]
                    ok = [r for r in every if r['status'] == 'ok']
                    if not every:
                        continue
                    several = [r for r in ok if r['k'] >= 2]
                    one = [r for r in ok if r['k'] == 1]
                    f.write('\t'.join([setting, str(n), method, str(len(every)), str(len(every) - len(ok))] + [
                        f'{mean([number(r[c]) for r in ok]):.4g}' for c in ('seconds', 'nodes', 'cycles', 'length', 'fit')] + [
                        f'{mean([r[c] for r in several]):.4f}' for c in ('recall', 'precision', 'conflictRecall', 'conflictPrecision', 'falseSplit', 'overDisplay')] + [
                        f'{mean([r["falseSplit"] for r in one]):.4f}']) + '\n')

    # numbers for the text: averages over all numbers of taxa
    macros = {}
    for setting, word in (('inf', 'Inf'), ('low', 'Low'), ('mid', 'Mid'), ('high', 'High')):
        for method, short in (('Distances', 'D'), ('NeighborNet', 'NN'), ('RazorNet', 'RN'), ('RazorNet+Simplify', 'RNS')):
            ok = [r for r in done if r['setting'] == setting and r['method'] == method]
            several = [r for r in ok if r['k'] >= 2]
            one = [r for r in ok if r['k'] == 1]
            name = 'simT' + word + short
            macros[name + 'Recall'] = f'{100 * mean([r["conflictRecall"] for r in several]):.0f}'
            macros[name + 'Precision'] = f'{100 * mean([r["conflictPrecision"] for r in several]):.0f}'
            macros[name + 'False'] = f'{100 * mean([r["falseSplit"] for r in several]):.1f}'
            macros[name + 'FalseOne'] = f'{100 * mean([r["falseSplit"] for r in one]):.1f}'
            macros[name + 'Over'] = f'{100 * mean([r["overDisplay"] for r in several]):.1f}'
            if method != 'Distances':
                macros[name + 'Cycles'] = f'{mean([number(r["cycles"]) for r in ok]):.0f}'
                macros[name + 'Length'] = f'{mean([number(r["length"]) for r in ok]):.0f}'
        timeouts = [r for r in rows if r['setting'] == setting and r['method'] == 'RazorNet+Simplify' and r['status'] != 'ok']
        macros['simT' + word + 'Timeouts'] = str(len(timeouts))
    counts = collections.Counter((r['setting'], r['n'], r['k']) for r in rows if r['method'] == 'NeighborNet')
    macros['simTReplicates'] = str(min(counts.values())) if counts else '0'
    macros['simTDatasets'] = str(sum(counts.values()))
    return macros


# ---------------------------------------------------------------------------------------------------------
# study 2: samples from a population
# ---------------------------------------------------------------------------------------------------------
METHODS_2 = ['MedianJoining', 'TCS', 'RazorHaploNet', 'RazorHaploNet+Simplify']
PANELS_2 = [('linkRecall', 'true links found'), ('linkPrecision', 'links that are true'),
            ('mutationRatio', 'mutations drawn\n(relative to the truth)'), ('exactPairs', 'pairs displayed exactly')]


def study2(rows, figure_directory, table_directory):
    for row in rows:
        row['n'], row['theta'] = (int(x) for x in row['cell'][1:].split('t'))
        row['linkRecall'] = ratio(row, 'sharedEdges', 'trueEdges')
        row['linkPrecision'] = ratio(row, 'sharedEdges', 'shownEdges')
        row['mutationRatio'] = ratio(row, 'mutations', 'trueMutations')
        row['exactPairs'] = number(row['exactPairs'])
        row['recall'] = ratio(row, 'sharedSupport', 'trueSupport')
        row['precision'] = ratio(row, 'sharedSupport', 'shownSupport')
    done = [r for r in rows if r['status'] == 'ok']
    sizes = sorted(set(r['n'] for r in rows))
    thetas = sorted(set(r['theta'] for r in rows))

    for n in sizes:
        figure, axes = plt.subplots(1, 4, figsize=(7.2, 2.5))
        for ax, (measure, label) in zip(axes, PANELS_2):
            for method in METHODS_2 + (['Truth'] if measure == 'exactPairs' else []):
                values = collections.defaultdict(list)
                for r in done:
                    if r['setting'] == 'finite' and r['n'] == n and r['method'] == method:
                        values[r['theta']].append(r[measure])
                if values:
                    xs = sorted(values)
                    draw_series(ax, xs, values, method, positions=[thetas.index(x) for x in xs], use_median=True)
            ax.set_title(label, fontweight='normal', fontsize=8)
            ax.set_xlabel(r'diversity $\theta$')
            ax.set_xticks(range(len(thetas)))
            ax.set_xticklabels([str(t) for t in thetas])
            style_axes(ax)
            if measure == 'mutationRatio':
                ax.set_ylim(bottom=0.9)
                ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:g}x'))
            else:
                ax.set_ylim(0, 1.05)
                ax.yaxis.set_major_formatter(FuncFormatter(percent))
        handles, labels = axes[3].get_legend_handles_labels()
        figure.legend(handles, labels, loc='upper center', ncol=5, frameon=False, bbox_to_anchor=(0.5, 1.0),
                      handlelength=2.2, columnspacing=1.6)
        figure.tight_layout(rect=(0, 0, 1, 0.9))
        figure.savefig(os.path.join(figure_directory, f'sim-haplotypes-n{n}.pdf'))
        plt.close(figure)

    with open(os.path.join(table_directory, 'simulation-study-2-summary.tsv'), 'w') as f:
        f.write('model\tsamples\ttheta\tmethod\tdatasets\ttimeouts\thaplotypes\tseconds\tauxiliary\tcycles\tmutations'
                '\ttrueMutations\texcess\texactPairs\tlinkRecall\tlinkPrecision\trecall\tprecision\n')
        for model in ('inf', 'finite'):
            for n in sizes:
                for theta in thetas:
                    for method in ['Truth'] + METHODS_2:
                        every = [r for r in rows if r['setting'] == model and r['n'] == n and r['theta'] == theta and r['method'] == method]
                        ok = [r for r in every if r['status'] == 'ok']
                        if not every:
                            continue
                        f.write('\t'.join([model, str(n), str(theta), method, str(len(every)), str(len(every) - len(ok))] + [
                            f'{median([number(r[c]) for r in ok]):.4g}' for c in ('haplotypes', 'seconds', 'auxiliary', 'cycles', 'mutations', 'trueMutations', 'excess')] + [
                            f'{median([r[c] for r in ok]):.4f}' for c in ('exactPairs', 'linkRecall', 'linkPrecision', 'recall', 'precision')]) + '\n')

    # numbers for the text: medians over all sample sizes (a few samples give very large networks)
    macros = {}
    words = {5: 'Five', 10: 'Ten', 20: 'Twenty', 50: 'Fifty'}
    for model, word in (('inf', 'Inf'), ('finite', 'Fin')):
        for theta in thetas:
            for method, short in (('Truth', 'True'), ('MedianJoining', 'MJ'), ('TCS', 'TCS'), ('RazorHaploNet', 'RHN'), ('RazorHaploNet+Simplify', 'RHNS')):
                ok = [r for r in done if r['setting'] == model and r['theta'] == theta and r['method'] == method]
                name = 'simH' + word + words[theta] + short
                macros[name + 'LinkRecall'] = f'{100 * median([r["linkRecall"] for r in ok]):.0f}'
                macros[name + 'LinkPrecision'] = f'{100 * median([r["linkPrecision"] for r in ok]):.0f}'
                macros[name + 'Mutations'] = f'{median([r["mutationRatio"] for r in ok]):.2f}'
                macros[name + 'Exact'] = f'{100 * median([r["exactPairs"] for r in ok]):.0f}'
                macros[name + 'Cycles'] = f'{median([number(r["cycles"]) for r in ok]):.0f}'
                macros[name + 'MaxCycles'] = f'{max([number(r["cycles"]) for r in ok], default=math.nan):.0f}'
        every = [r for r in rows if r['setting'] == model and r['method'] == 'MedianJoining']
        macros['simH' + word + 'MJTimeouts'] = str(len([r for r in every if r['status'] != 'ok']))
        # samples for which a method does not return the true haplotype tree
        for method, short in (('MedianJoining', 'MJ'), ('TCS', 'TCS'), ('RazorHaploNet', 'RHN'), ('RazorHaploNet+Simplify', 'RHNS')):
            ok = [r for r in done if r['setting'] == model and r['method'] == method]
            macros['simH' + word + short + 'NotTheTree'] = str(len([r for r in ok if r['linkRecall'] < 1 or r['linkPrecision'] < 1]))
            macros['simH' + word + short + 'Samples'] = str(len(ok))
    counts = collections.Counter((r['setting'], r['n'], r['theta']) for r in rows if r['method'] == 'TCS')
    macros['simHReplicates'] = str(min(counts.values())) if counts else '0'
    macros['simHDatasets'] = str(sum(counts.values()))
    return macros


def main():
    study1_table, study2_table, figure_directory, table_directory = sys.argv[1:5]
    os.makedirs(figure_directory, exist_ok=True)
    os.makedirs(table_directory, exist_ok=True)
    macros = study1(read(study1_table), figure_directory, table_directory)
    macros.update(study2(read(study2_table), figure_directory, table_directory))
    # the numbers that the text of the manuscript quotes, as macros
    with open(os.path.join(table_directory, 'simulation-numbers.tex'), 'w') as f:
        f.write('% written by experiments/figures.py: means (study 1) and medians (study 2) over the simulated data sets\n')
        for name in sorted(macros):
            f.write(f'\\newcommand{{\\{name}}}{{{macros[name]}}}\n')


if __name__ == '__main__':
    main()
