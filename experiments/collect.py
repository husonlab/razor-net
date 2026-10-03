#!/usr/bin/env python3
"""Collects the per-data-set reports of a simulation study into one table.

    collect.py <directory of the study> <output.tsv>

The study directory holds <setting>/<cell>/rep<NNN>/report.tsv, as written by run-simulation-study-*.sh. The
output has the rows of all reports, with the setting, the cell and the replicate as the first three columns and
the path of the input removed.
"""
import csv
import glob
import os
import sys


def main():
    root, output = sys.argv[1], sys.argv[2]
    rows, header = [], None
    for path in sorted(glob.glob(os.path.join(root, '*', '*', 'rep*', 'report.tsv'))):
        setting, cell, replicate = path.split(os.sep)[-4:-1]
        with open(path, newline='') as f:
            for row in csv.DictReader(f, delimiter='\t'):
                if header is None:
                    header = ['setting', 'cell', 'replicate'] + [c for c in row if c != 'file']
                row.update(setting=setting, cell=cell, replicate=replicate)
                rows.append(row)
    with open(output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=header, delimiter='\t', extrasaction='ignore', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    print(f'{len(rows)} rows from {len(set((r["setting"], r["cell"], r["replicate"]) for r in rows))} data sets -> {output}')


if __name__ == '__main__':
    main()
