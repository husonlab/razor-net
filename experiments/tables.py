#!/usr/bin/env python3
r"""Writes the LaTeX rows of the supplementary tables from the results.

    tables.py <results directory> <table directory>

Reads published-alignments.tsv, simulation-study-1-summary.tsv and simulation-study-2-summary.tsv, and writes
alignments-rows.tex, sim-trees-rows.tex and sim-haplotypes-rows.tex: the rows of a table, ending with its bottom
rule, to be read by \input inside the tabular.
"""
import csv
import os
import sys

METHOD = {'NeighborNet': r'\NN', 'RazorNet': r'\RN', 'RazorNet+Simplify': r'\RN, Simplify', 'Distances': 'distances',
          'MedianJoining': r'\MJ', 'TCS': r'\TCS', 'RazorHaploNet': r'\RHN', 'RazorHaploNet+Simplify': r'\RHN, Simplify',
          'Truth': 'true tree'}
SETTING = {'inf': 'infinite sites', 'low': 'low', 'mid': 'medium', 'high': 'high', 'finite': 'finite sites'}


def read(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def tex(text):
    return text.replace('_', r'\_')


def percent(text):
    return '--' if text in ('', 'nan') else f'{100 * float(text):.1f}'


def plain(text, digits=0):
    return '--' if text in ('', 'nan') else f'{float(text):.{digits}f}'


def alignments(results, tables):
    by_file = {}
    for row in read(os.path.join(results, 'published-alignments.tsv')):
        by_file.setdefault(row['file'], {})[row['method']] = row
    lines = []
    sources = ('rhabdoviruses', 'XuEtAl2026', 'SalugaEtAl2022', 'Issakka_ITSClustering')  # the order of the rows
    for file, methods in sorted(by_file.items(), key=lambda item: [source in item[0] for source in sources].index(True)):
        name = os.path.basename(file)
        for old in ('ALIGNMENT_', 'ALIGNEMNT_', '_Nucleotide alignment', '.fasta', '.fas', '.phy', '_SEQ'):
            name = name.replace(old, '')
        first = next(iter(methods.values()))
        cells = []
        for method in ('MedianJoining', 'TCS', 'RazorHaploNet', 'RazorHaploNet+Simplify'):
            row = methods[method]
            cells.append(f"{row['mutations']} & {row['cycles']} & {row['excess']}")
        lines.append(rf"\texttt{{{tex(name)}}} & {first['sequences']} & {first['haplotypes']} & {first['sites']} & "
                     + ' & '.join(cells) + r' \\')
    with open(os.path.join(tables, 'alignments-rows.tex'), 'w') as f:
        f.write('\n'.join(lines) + '\n\\bottomrule\n')


def sim_trees(results, tables):
    lines, previous = [], None
    for row in read(os.path.join(results, 'simulation-study-1-summary.tsv')):
        key = (row['setting'], row['taxa'])
        if previous is not None and key != previous:
            lines.append(r'\addlinespace')
        label = f"{SETTING[row['setting']]} & {row['taxa']}" if key != previous else ' & '
        previous = key
        size = (f"{plain(row['nodes'])} & {plain(row['cycles'])} & {plain(row['length'])} & {plain(row['fit'], 2)}"
                if row['method'] != 'Distances' else ' & & & ')
        lines.append(f"{label} & {METHOD[row['method']]} & {size} & {percent(row['conflictRecall'])} & "
                     f"{percent(row['conflictPrecision'])} & {percent(row['falseSplit'])} & "
                     f"{percent(row['falseSplitOneTree'])} & {row['timeouts']} \\\\")
    with open(os.path.join(tables, 'sim-trees-rows.tex'), 'w') as f:
        f.write('\n'.join(lines) + '\n\\bottomrule\n')


def sim_haplotypes(results, tables):
    lines, previous = [], None
    for row in read(os.path.join(results, 'simulation-study-2-summary.tsv')):
        if row['model'] != 'finite':
            continue
        key = (row['samples'], row['theta'])
        if previous is not None and key != previous:
            lines.append(r'\addlinespace')
        label = f"{row['samples']} & {row['theta']} & {plain(row['haplotypes'], 1)}" if key != previous else ' & & '
        previous = key
        lines.append(f"{label} & {METHOD[row['method']]} & {plain(row['auxiliary'], 1)} & {plain(row['cycles'], 1)} & "
                     f"{plain(row['mutations'], 1)} & {plain(row['excess'])} & {percent(row['exactPairs'])} & "
                     f"{percent(row['linkRecall'])} & {percent(row['linkPrecision'])} & {row['timeouts']} \\\\")
    with open(os.path.join(tables, 'sim-haplotypes-rows.tex'), 'w') as f:
        f.write('\n'.join(lines) + '\n\\bottomrule\n')


def main():
    results, tables = sys.argv[1:3]
    os.makedirs(tables, exist_ok=True)
    alignments(results, tables)
    sim_trees(results, tables)
    sim_haplotypes(results, tables)


if __name__ == '__main__':
    main()
