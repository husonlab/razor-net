# Checking our TCS construction against PopART

The seven alignments of Table 1 of the manuscript, each as the common input of the comparison (the alignment
without every column that contains a gap, a missing or an ambiguous character; identical sequences are left
for PopART to merge), in the NEXUS dialect that PopART reads. Written by

```bash
java -cp "$JAVA_CP" razornet.testsuite.PopartCheck --out examples/popart-check <alignments>
```

which also prints the size of the TCS network that our construction (`razornet.algorithms.tcs.TCSConstruction`,
written after PopART's `TCS.cpp`) builds on the same input. Because the columns are already removed, PopART's
own masking has nothing to do, and the network it draws should have the same size as ours.

## What to do

Open each file in PopART, choose **Network → TCS Network**, and compare with the table below:

- **intermediates**: the small unlabeled nodes of the drawing, inferred intermediates where three or more
  edges meet;
- **steps**: the total number of mutations drawn, the hatch marks summed over all edges (a direct edge between
  two haplotypes without a hatch mark counts one);
- **cycles**: the number of independent cycles, edges minus nodes plus one.

The numbers are those of the construction, before sequences are assigned to the intermediate nodes. Table 1
reports the *labeled* TCS network, which can differ: the labeling gives every intermediate a sequence and every
edge its mutations, and then removes intermediates and edges that have become redundant, so the "nodes added"
and "mutations" of Table 1 are not the numbers to compare with PopART. They are given in the last columns for
reference.

| data set | sequences | haplotypes | sites | intermediates | steps | cycles | Table 1, labeled: nodes added / cycles / mutations |
|---|---:|---:|---:|---:|---:|---:|---|
| CCapSV | 19 | 15 | 1277 | 11 | 52 | 3 | 9 / 3 / 52 |
| DImmSV | 87 | 43 | 1277 | 3 | 70 | 1 | 3 / 1 / 70 |
| PAegRV | 131 | 90 | 1246 | 39 | 401 | 9 | 31 / 8 / 423 |
| Diamond Princess | 482 | 257 | 444 | 27 | 489 | 12 | 24 / 12 / 495 |
| *Syntrichia*, ITS, haplotype network (116) | 116 | 16 | 712 | 5 | 65 | 3 | 4 / 3 / 66 |
| *Syntrichia*, ITS, phylogeography (127) | 127 | 13 | 410 | 11 | 177 | 3 | 9 / 3 / 190 |
| *Syntrichia*, maritime Antarctic (49) | 49 | 8 | 1757 | 1 | 13 | 0 | 1 / 0 / 13 |

For the 116-sequence ITS alignment the construction's numbers are those of Figure 1a of the manuscript (65
mutations, 3 cycles, 5 intermediates), which is the topology the authors published.

## Result of the check (DH, 2026-10-24)

PopART's TCS networks on the seven `…-popart.nex` files have the numbers of intermediates and cycles of the
table above on all seven alignments, and the numbers of steps on the four whose steps could be counted
(CCapSV 52, DImmSV 70, *Syntrichia* ITS 116 sequences 65, maritime Antarctic 13); the steps of PAegRV (401),
Diamond Princess (489) and *Syntrichia* ITS 127 sequences (177) were not counted. Median joining in PopART and
in RazorNet agreed on the same inputs as well. On the two ITS alignments a network computed in RazorNet's
window differed, for the reason given in the next section.

## The binary indel codes of the two *Syntrichia* ITS alignments

The 116-sequence and the 127-sequence ITS alignments carry, after the nucleotides, 13 and 21 columns of 0 and
1 in which the authors coded indels. Three programs read them differently:

- **PopART counts them as characters.** For DNA data its `Sequence::isAmbiguousChar` treats only the gap and
  the IUPAC codes as ambiguous; any other symbol is a state, and a difference in a 0/1 column is a step
  (`HapNet.cpp`, `pairwiseDistance`).
- **Our comparison counts them as states too** (`CommonInput`: a gap, a missing character or an IUPAC code
  makes a column incomplete, every other symbol is a state), so the expected numbers above include them, as do
  the tables of the paper.
- **RazorNet's own import of these files ignores them.** The CHARACTERS block says `DATATYPE=DNA`, and the
  NEXUS reader of SplitsTree drops the sites with 0 and 1 with a warning. A network computed in RazorNet on one
  of these files is therefore on 699 or 389 sites, not 712 or 410, and is not the network to compare with
  PopART. On the 127-sequence alignment it gives, for TCS, 162 mutations, an excess of 125 and 4 cycles
  (labeled), where the paper's input gives 190, 140 and 3.

For both readings, the construction's numbers are:

| data set | sites | intermediates | steps | cycles |
|---|---:|---:|---:|---:|
| *Syntrichia* ITS, 116 sequences, with the 13 indel columns | 712 | 5 | 65 | 3 |
| the same, nucleotides only (`…-dna-only-popart.nex`): 15 haplotypes, two merge | 699 | 2 | 37 | 4 |
| *Syntrichia* ITS, 127 sequences, with the 21 indel columns | 410 | 11 | 177 | 3 |
| the same, nucleotides only (`…-dna-only-popart.nex`): 13 haplotypes | 389 | 8 | 158 | 4 |

PopART on `…-popart.nex` is expected to match the first row of each pair. The `…-dna-only-popart.nex` files
hold the nucleotide columns alone; on them, PopART and a network computed in RazorNet see the same input.

| file | alignment |
|---|---|
| `CCapSV-popart.nex`, `DImmSV-popart.nex`, `PAegRV-popart.nex` | `examples/rhabdoviruses/*.fas` (Longdon et al. 2017) |
| `Diamond_Princess-popart.nex` | `examples/publications/XuEtAl2026/Diamond_Princess.phy` (Xu et al. 2026) |
| `ALIGNMENT_ITS_Haplotype_Network_116_SEQ-popart.nex` | Saługa et al. 2022, Dryad deposit |
| `ALIGNMENT_ITS_Phylogeography_127_SEQ-popart.nex` | Saługa et al. 2022, Dryad deposit |
| `…-dna-only-popart.nex` (two files) | the two ITS alignments without their 0/1 indel columns, see above |
| `ALIGNMENT_Maritime_Antarctic_Population_Genetics_49_SEQ-popart.nex` | Saługa et al. 2022, Dryad deposit |

The two files `…-popart.nex` elsewhere in `examples/publications/` (the 116 ITS sequences and the Diamond
Princess genomes, with traits) hold the full alignments for drawing; PopART masks their columns itself.
