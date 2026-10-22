#!/usr/bin/env python3
"""
Write the ITS alignment of Saługa et al. (2022) as a NEXUS file for PopART.

Reads  doi_10_5061_dryad_h18931zp3__v20220812/ALIGNMENT_ITS_Haplotype_Network_116_SEQ.fasta
       (the matrix the authors ran in PopART 1.7, 116 sequences, 951 sites, as deposited)
and    ALIGNMENT_ITS_Haplotype_Network_116_SEQ.stree6, whose TAXINFO lists the members of
       each ITS haplotype of the published Figure 2b (R02 ... R18; PAT-05 within R17, as
       PopART groups it). The three haplotypes of a single sequence (R04, R10, R11) name no
       member there; they are matched to the alignment by their sequence in the file's
       CHARACTERS block.

Writes ALIGNMENT_ITS_Haplotype_Network_116_SEQ-popart.nex: TAXA, CHARACTERS and PopART's
TRAITS block, one trait per haplotype, so that PopART colours each node as the published
figure does. In the shape of PopART's example file: no comments, nothing quoted, one space
between a name and its data, names of letters, digits and '_' only ('-' becomes '_').
Daniel Huson / AI agent, 10.2026
"""

import os
import re
import sys
from collections import OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
FASTA = os.path.join(HERE, "doi_10_5061_dryad_h18931zp3__v20220812", "ALIGNMENT_ITS_Haplotype_Network_116_SEQ.fasta")
WORKFLOW = os.path.join(HERE, "ALIGNMENT_ITS_Haplotype_Network_116_SEQ.stree6")
OUT = os.path.join(HERE, "ALIGNMENT_ITS_Haplotype_Network_116_SEQ-popart.nex")


def read_fasta(path):
    seqs = OrderedDict()
    name = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].strip()
                seqs[name] = ""
            else:
                seqs[name] += line
    lengths = {len(s) for s in seqs.values()}
    if len(lengths) != 1:
        sys.exit("sequences of unequal length: %s" % sorted(lengths))
    return seqs, lengths.pop()


def read_haplotypes(path):
    """haplotype label -> members, from the TAXLABELS and TAXINFO of the first TAXA block"""
    text = open(path).read()
    block = text[text.index("BEGIN TAXA;"):text.index("END; [TAXA]")]
    labels = re.findall(r"\[\d+\]\s+'([^']*)'", block[block.index("TAXLABELS"):block.index("TAXINFO")])
    infos = re.findall(r"\[\d+\]\s+'([^']*)'", block[block.index("TAXINFO"):])
    if len(labels) != len(infos):
        sys.exit("TAXLABELS and TAXINFO differ in length")
    return OrderedDict((label, info.split(",")) for label, info in zip(labels, infos))


def read_workflow_sequences(path):
    """label -> sequence, from the first CHARACTERS block of the workflow file (interleaved)"""
    text = open(path).read()
    start = text.index("BEGIN CHARACTERS;")
    block = text[start:text.index("END; [CHARACTERS]", start)]
    rows = OrderedDict()
    for line in block[block.index("MATRIX") + len("MATRIX"):].split("\n"):
        line = line.strip()
        if not line or line == ";":
            continue
        label, chunk = line.split(None, 1)
        rows[label] = rows.get(label, "") + chunk.strip()
    return rows


def popart_label(name):
    return re.sub(r"[^A-Za-z0-9_]", "_", name)


def main():
    seqs, nchar = read_fasta(FASTA)
    haplotypes = read_haplotypes(WORKFLOW)
    workflow_seqs = read_workflow_sequences(WORKFLOW)
    haplotype_of = {}
    for label, members in haplotypes.items():
        if members == [label]:   # a haplotype of one sequence: find that sequence
            members = [n for n, q in seqs.items() if q.upper() == workflow_seqs[label].upper()]
            if len(members) != 1:
                sys.exit("%s matches %d sequences of the alignment" % (label, len(members)))
            haplotypes[label] = members
        for member in members:
            if member in haplotype_of:
                sys.exit("%s listed under %s and %s" % (member, haplotype_of[member], label))
            haplotype_of[member] = label
    missing = [n for n in seqs if n not in haplotype_of]
    extra = [m for m in haplotype_of if m not in seqs]
    if missing or extra:
        sys.exit("sequences without a haplotype: %s; members not in the alignment: %s" % (missing, extra))
    traits = sorted(haplotypes)
    names = list(seqs)
    popart_names = [popart_label(n) for n in names]
    if len(set(popart_names)) != len(popart_names):
        sys.exit("two sequence names coincide after replacing characters by '_'")

    with open(OUT, "w") as w:
        w.write("#NEXUS\n\n")
        w.write("BEGIN TAXA;\n")
        w.write("DIMENSIONS NTAX=%d;\n" % len(names))
        w.write("TAXLABELS\n")
        for pname in popart_names:
            w.write("%s\n" % pname)
        w.write(";\nEND;\n\n")
        w.write("BEGIN CHARACTERS;\n")
        w.write("DIMENSIONS NCHAR=%d;\n" % nchar)
        w.write("FORMAT DATATYPE=DNA MISSING=? GAP=-;\n")
        w.write("MATRIX\n")
        for name, pname in zip(names, popart_names):
            w.write("%s %s\n" % (pname, seqs[name].upper()))
        w.write(";\nEND;\n\n")
        w.write("BEGIN TRAITS;\n")
        w.write("Dimensions NTRAITS=%d;\n" % len(traits))
        w.write("Format labels=yes missing=? separator=Comma;\n")
        w.write("TraitLabels " + " ".join(traits) + ";\n")
        w.write("Matrix\n")
        for name, pname in zip(names, popart_names):
            w.write("%s %s\n" % (pname, ",".join("1" if haplotype_of[name] == t else "0" for t in traits)))
        w.write(";\nEND;\n")

    alphabet = sorted(set("".join(seqs.values()).upper()))
    print("sequences %d, sites %d, alphabet %s" % (len(names), nchar, "".join(alphabet)))
    print("haplotypes: " + ", ".join("%s (%d)" % (t, len(haplotypes[t])) for t in traits))
    print("single sequences: " + ", ".join("%s = %s" % (t, haplotypes[t][0]) for t in traits if len(haplotypes[t]) == 1))
    print("wrote " + OUT)


if __name__ == "__main__":
    main()
