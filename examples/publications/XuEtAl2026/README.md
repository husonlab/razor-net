# Xu et al. (2026) — Diamond Princess SARS-CoV-2 haplotypes

SARS-CoV-2 sequences from the Diamond Princess cruise-ship outbreak together with
contemporaneous sequences from elsewhere, as distributed with:

> Bo Xu, Lun Li, Cuiping Li, Anke Wang, Zhuojing Fan, Shuhui Song. *HaploThread: a
> scalable desktop tool for efficient haplotype network inference and interactive
> visualization.* Molecular Biology and Evolution **43**(3), March 2026, msag052.
> <https://doi.org/10.1093/molbev/msag052>

482 sequences over 444 variable sites (alphabet `{A,C,G,T}` — no gaps, no ambiguity
codes), collected 2019-12-24 … 2020-02-19. Collapsing identical sequences gives **257
unique haplotypes**.

## Files

| File | What it is |
|---|---|
| `Diamond_Princess.phy` | *as received* — relaxed PHYLIP, 482 × 444 |
| `Diamond_Princess.tsv` | *as received* — `*`, name, date, location, `*`, `*`; no header |
| `Diamond_Princess.vcf` | *as received* — the same variants as VCF (not used below) |
| **`Diamond_Princess.nex`** | **self-contained NEXUS: TAXA + CHARACTERS + TRAITS. Start here.** |
| `Diamond_Princess-traits-matrix.tsv` | the traits as a matrix with a `Coordinates` row, for routes 2 and 3 below |
| `Diamond_Princess-popart.nex` | the same sequences and the eight location traits in the NEXUS dialect that PopART reads: no `ntax` in the CHARACTERS block (PopART refuses it: "ntax only allowed in characters block along with new taxa"), PopART's own TRAITS block in the shape of PopART's example file (nothing quoted, no comments, one space between a name and its data, no trait coordinates), `/`, `-` and `.` in sequence names replaced by `_`; the month bins are left out because PopART draws every trait of the block into one pie |
| `make-traits.py` | regenerates the three files above from the `.phy` and `.tsv` |

### Workflow files

| File | What it holds |
|---|---|
| `Diamond_Princess-four-razornet.stree6` | the 257 haplotypes with their traits; RazorHaploNet (default options), then the Shave filter at its default target, Balanced: the network of Figure 2 of the manuscript, drawn with the sampling locations as a legend (the months are not shown) |
| `Diamond_Princess-median-joining.stree6` | the same haplotypes under median joining |
| `Diamond_Princess-four.stree6` | RazorHaploNet (100 restarts), median joining, TCS, and RazorHaploNet with the Shave filter, side by side |
| `Diamond_Princess.stree6` | RazorHaploNet (100 restarts), median joining, and RazorHaploNet with the Stretch filter |

## Three ways to load it

1. **Open `Diamond_Princess.nex`.** Everything in one file, nothing to collapse or
   re-import. Keeps the exact per-sample `Date` (see below).
2. **File → Open Haplotype Data**, sequences = `Diamond_Princess.phy`, traits =
   `Diamond_Princess-traits-matrix.tsv`. This collapses the 482 sequences to the 257
   unique haplotypes and pools the trait counts, then runs the chosen network method.
   Needs a SplitsTree built after 2026-08-20 — see below.
3. **Open `Diamond_Princess.phy`, then File → Import → Taxon Traits** with the same
   `Diamond_Princess-traits-matrix.tsv`, if you would rather start from the raw
   alignment. This does not collapse anything, so you keep all 482 taxa.

## Traits

`TraitsBlock` is trait-*per-column*, and `TRAITLATITUDE` / `TRAITLONGITUDE` give one
coordinate per **trait**, not per taxon. Each location is therefore its own trait, and
the matrix holds counts.

| # | Trait | Kind | Latitude | Longitude | n |
|---|---|---|---:|---:|---:|
| 1 | `China` | numeric | 35.8617 | 104.1954 | 244 |
| 2 | `Japan` | numeric | 36.2048 | 138.2529 | 9 |
| 3 | `DP` | numeric | 35.4500 | 139.6667 | 72 |
| 4 | `CruiseA` | numeric | 38.2627 | −121.9272 | 12 |
| 5 | `other_Asia` | numeric | 20.9259 | 111.7914 | 100 |
| 6 | `North_America` | numeric | 42.5453 | −99.8740 | 18 |
| 7 | `Europe` | numeric | 47.5108 | 4.9789 | 16 |
| 8 | `Oceania` | numeric | −25.2744 | 133.7751 | 11 |
| 9 | `2019-12` | numeric | — | — | 17 |
| 10 | `2020-01` | numeric | — | — | 252 |
| 11 | `2020-02` | numeric | — | — | 213 |

Every sample contributes exactly 1 to the location group and 1 to the month group.
`make-traits.py --with-date` adds a 12th, categorical `Date` trait — see below.

### Where the coordinates come from

- **`DP`** — the Diamond Princess's quarantine berth, **Port of Yokohama** (Daikoku
  Pier), Japan. Note that `Japan/Hu_DP_Kng_19-031/2020` carries `DP` in its name but is
  labelled `Japan` in the source `.tsv`; the source labelling is left untouched.
- **`CruiseA`** — the 12 `USA/CruiseA-*` samples are GenBank MT159705ff., submitted by
  the CDC under the title *"Cruise A sequences"*, `geo_loc_name="USA"`, collected
  2020-02-17/18: US citizens repatriated from the Diamond Princess on the 16–17 February
  charter flights and sampled after arrival. Placed at **Travis AFB, California**, the
  main repatriation destination, which keeps them visually distinct from `DP`. Move them
  to Yokohama (35.4500, 139.6667) if you would rather show where the infections were
  acquired.
- **`China`** and **`Japan`** — country centroids. The Chinese samples are genuinely
  country-wide (Shanghai 42, Hangzhou 35, Guangdong 31, Wuhan 26, Guangzhou 21,
  Beijing 20, Jiangxi, Shandong, …), so a Wuhan point would misrepresent them.
- **`other_Asia`**, **`North_America`**, **`Europe`**, **`Oceania`** — the
  **sample-weighted centroid** of their member countries, rather than a hand-picked
  regional point, so each pie sits at its own centre of mass:

  | Region | Composition |
  |---|---|
  | `other_Asia` | HongKong ×29, Korea ×14, Thailand ×12, Singapore ×11, SouthKorea ×9, Taiwan ×8, Malaysia ×6, Vietnam ×5, India ×2, UnitedArabEmirates ×2, Nepal ×1, Cambodia ×1 |
  | `North_America` | USA ×15, Canada ×3 |
  | `Europe` | France ×10, Italy ×2, Germany ×1, England ×1, Belgium ×1, Sweden ×1 |
  | `Oceania` | Australia ×11 |

  `other_Asia` therefore lands just off the Guangdong coast — the honest centre of a
  group dominated by Hong Kong, Korea and Southeast Asia.

To change any of this, edit `POINT` / `COUNTRY` in `make-traits.py` and re-run it.

### Date

- **`2019-12` / `2020-01` / `2020-02`** — *numeric* month bins, shipped by default.
  They survive haplotype collapsing, so they give a temporal read-out on every route.
- **`Date`** — an optional *categorical* trait holding the exact value (`2020-01-25`, …;
  48 distinct values, one of which is the partial `2020-01`, as in the source). Add it
  with `python3 make-traits.py --with-date`. Categorical traits are ignored by the pie
  charts and by the World Map, so the column costs nothing visually; use
  **Select Traits → Date → value** to pick out a day. It is *dropped* when identical
  haplotypes are collapsed, so it never reaches the Open-Haplotype-Data route.

  It is **off by default** because of the third point below: on a SplitsTree built before
  2026-08-20 this one column disables the Taxon Traits Pie Charts for the whole dataset.
  Turn it on once you are running a newer build.

## Three things to know

- **The pies mix both groups by default.** With all traits active, each node's pie sums
  to 2 — half location, half month. Deselect one group under **Traits → Active** to get
  a clean location pie or a clean month pie. The World Map is unaffected: it draws only
  traits that have coordinates, so it always shows locations alone.
- **You need a SplitsTree built after 2026-08-20 for the World Map.** Preparing this
  dataset turned up an off-by-one in `WorldMapView`: `updateNodes` and `computeMaxCount`
  both looped `traitId = 1; traitId < traitsBlock.getNTraits()` while traits are 1-based,
  so the last trait was never drawn and never counted towards the pie scale — with *n*
  locations you silently saw *n*−1 pies. Fixed on 2026-08-20 (`WorldMapView.java:321`
  and `:362`). All three files here are ordered so that the last trait is a
  coordinate-free one anyway, so they were never affected and still work on an older
  build.
- **Route 2 needs a SplitsTree built after 2026-08-20.** The Open Haplotype Data dialog
  used to read only a one-line-per-taxon traits layout, in which a line's single
  `(lat,long)` applies to *every* trait named on that line — so the month bins could not
  be left unplaced, and a separate file in that layout had to be shipped alongside this
  one. Handed the matrix file, the old dialog skipped the `traits` and `Coordinates`
  rows, read each matrix *cell* as a trait name, and silently produced two traits called
  `0` and `1` with no coordinates and no world map. It now recognizes the matrix layout
  by its `traits` keyword and shares one parser with route 3, so a single traits file
  serves both.
- **A categorical trait used to disable every pie chart.**
  `TraitsBlock.setInducedTraits` — the copy the Taxa Filter makes for the working taxa —
  copied labels with `getTraitValueLabel`, which never returns null and formats the
  number instead. So one categorical trait anywhere in the block turned *every* trait of
  the working copy categorical, `getNumberNumericalTraits()` fell to 0, and
  `TraitsFormatPresenter` greyed out the **Taxon Traits Pie Charts** pane. Fixed on
  2026-08-20 (`TraitsBlock.java:326-334`); until you are on a build with that fix, leave
  `--with-date` off.
