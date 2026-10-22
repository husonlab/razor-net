# Saługa et al. (2022) — ITS and cpDNA haplotype networks of *Syntrichia*

Marta Saługa, Ryszard Ochyra, Michał Ronikier. *Phylogeographical break and limited connectivity between
multiple refugia in panantarctic moss species.* Journal of Biogeography 49 (2022) 1991–2004,
<https://doi.org/10.1111/jbi.14476>. The paper is `2022_JBiogeogr_Saluga-etal_Syntrichia.pdf`; the data are
the Dryad deposit `doi_10_5061_dryad_h18931zp3__v20220812/` (Dryad 10.5061/dryad.h18931zp3, with its own
README; the `.hamming.dist` files next to the alignments are ours, written by `CharactersToDistances`).

Figure 1 of the manuscript is the ITS haplotype network of 116 sequences (`ALIGNMENT_ITS_Haplotype_Network_116_SEQ`,
the matrix the authors ran in PopART 1.7, Figure 2b of the paper).

| file | what it is |
|---|---|
| `ALIGNMENT_ITS_Haplotype_Network_116_SEQ.stree6` | the 116 sequences grouped into the 16 ITS haplotypes of the paper (R02–R18; its TAXINFO lists the members), with the traits and display labels below, RazorHaploNet and TCS |
| `…-Figure1a-PopART-TCS.stree6`, `…-haplotypes-fig1b.stree6`, `…-haplotypes-fig1c.stree6` | the three panels of Figure 1: TCS as PopART draws it, RazorHaploNet, RazorHaploNet with the Shave filter |
| `ALIGNMENT_cpDNA_Haplotype_Network_111_SEQ.stree6` | the cpDNA network of 111 sequences (Figure 2 of the paper, part c) |
| `traits.txt`, `display-labels.txt` | one trait per haplotype and the display labels, for File → Import in SplitsTree |
| `ALIGNMENT_ITS_Haplotype_Network_116_SEQ-popart.nex` | the 116 sequences as deposited, with one trait per haplotype, in the NEXUS dialect PopART reads; written by `make-popart.py` from the Dryad alignment and the TAXINFO of the workflow file (sequence names with `-` replaced by `_`) |
