#!/usr/bin/env python3
"""
Build the SplitsTree trait files for the Xu et al. (2026) Diamond Princess dataset.

Reads  Diamond_Princess.phy  (relaxed PHYLIP, 482 x 444, DNA)
and    Diamond_Princess.tsv  (no header; columns: *, name, date, location, *, *)

Writes Diamond_Princess.nex                  TAXA + CHARACTERS + TRAITS, self-contained
       Diamond_Princess-traits-matrix.tsv    trait matrix with a Coordinates row; read by
                                             the Open Haplotype Data dialog and by
                                             File -> Import -> Taxon Traits alike
       Diamond_Princess-popart.nex           the same sequences and the eight location
                                             traits in the NEXUS dialect PopART reads

By default the exact collection date is carried only by the three numerical month-bin
traits. Pass --with-date to add the categorical 'Date' trait holding the exact value per
sample; that needs SplitsTree built after 2026-08-20, because TraitsBlock.setInducedTraits
used to turn *every* trait of a block categorical as soon as one of them was, which left
the Taxon Traits Pie Charts pane with nothing numerical to draw and so disabled.

To move a location, edit POINT (or the COUNTRY table a region is averaged from) and re-run.
Daniel Huson / AI agent, 8.2026
"""

import os
import re
import sys
from collections import Counter, OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
PHY = os.path.join(HERE, "Diamond_Princess.phy")
TSV = os.path.join(HERE, "Diamond_Princess.tsv")

# --- coordinates -----------------------------------------------------------------
#
# Point locations.  DP is the Diamond Princess's quarantine berth (Daikoku Pier,
# Port of Yokohama).  CruiseA is GenBank MT159705ff. ("Cruise A sequences", CDC,
# geo_loc_name=USA, collected 2020-02-17/18) -- US repatriates flown off the Diamond
# Princess on the 16-17 Feb charters; placed at their repatriation destination.
POINT = {
    "China":   (35.8617, 104.1954),   # country centroid; samples span Wuhan, Shanghai,
                                      # Hangzhou, Guangdong, Beijing, Jiangxi, ...
    "Japan":   (36.2048, 138.2529),   # country centroid (Aichi, Tokyo, Kyoto, Osaka, Kanagawa)
    "DP":      (35.4500, 139.6667),   # Port of Yokohama, Japan
    "CruiseA": (38.2627, -121.9272),  # Travis AFB, California, USA
}

# Regional locations get a sample-weighted centroid of their member countries, keyed
# by the country prefix of the strain name.
COUNTRY = {
    "HongKong": (22.3193, 114.1694),
    "Korea": (35.9078, 127.7669),
    "SouthKorea": (35.9078, 127.7669),
    "Thailand": (15.8700, 100.9925),
    "Singapore": (1.3521, 103.8198),
    "Taiwan": (23.6978, 120.9605),
    "Malaysia": (4.2105, 101.9758),
    "Vietnam": (14.0583, 108.2772),
    "UnitedArabEmirates": (23.4241, 53.8478),
    "India": (20.5937, 78.9629),
    "Nepal": (28.3949, 84.1240),
    "Cambodia": (12.5657, 104.9910),
    "USA": (39.8283, -98.5795),
    "Canada": (56.1304, -106.3468),
    "France": (46.2276, 2.2137),
    "Italy": (41.8719, 12.5674),
    "Sweden": (60.1282, 18.6435),
    "Germany": (51.1657, 10.4515),
    "England": (52.3555, -1.1743),
    "Belgium": (50.5039, 4.4699),
    "Australia": (-25.2744, 133.7751),
}

REGIONS = ("other_Asia", "North_America", "Europe", "Oceania")

LOCATION_ORDER = ("China", "Japan", "DP", "CruiseA",
                  "other_Asia", "North_America", "Europe", "Oceania")


def country_of(name):
    """The country/city prefix of a strain name, e.g. 'France/IDF0372/2020' -> 'France'."""
    return name.split("/")[0]


def read_phylip(path):
    with open(path) as fh:
        header = fh.readline().split()
        ntax, nchar = int(header[0]), int(header[1])
        seqs = OrderedDict()
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            name, seq = line.split("\t", 1)
            seq = seq.strip()
            if len(seq) != nchar:
                sys.exit("length %d != %d for %s" % (len(seq), nchar, name))
            seqs[name] = seq
    if len(seqs) != ntax:
        sys.exit("read %d sequences, header says %d" % (len(seqs), ntax))
    return ntax, nchar, seqs


def read_traits(path):
    rows = OrderedDict()
    with open(path) as fh:
        for line in fh:
            tok = line.rstrip("\n").split("\t")
            if len(tok) < 4:
                continue
            rows[tok[1]] = (tok[2], tok[3])   # name -> (date, location)
    return rows


def weighted_centroid(names):
    """Sample-weighted mean of the member-country coordinates. All longitudes here
    lie well away from the antimeridian, so a plain mean is safe."""
    counts = Counter(country_of(n) for n in names)
    missing = [c for c in counts if c not in COUNTRY]
    if missing:
        sys.exit("no coordinates for: %s" % ", ".join(sorted(missing)))
    total = sum(counts.values())
    lat = sum(COUNTRY[c][0] * n for c, n in counts.items()) / total
    lon = sum(COUNTRY[c][1] * n for c, n in counts.items()) / total
    return round(lat, 4), round(lon, 4), counts


def main():
    ntax, nchar, seqs = read_phylip(PHY)
    meta = read_traits(TSV)

    names = list(seqs.keys())
    if set(names) != set(meta.keys()):
        sys.exit("taxon labels of .phy and .tsv differ")

    # coordinates per location
    coords = dict(POINT)
    breakdown = {}
    for region in REGIONS:
        members = [n for n in names if meta[n][1] == region]
        lat, lon, counts = weighted_centroid(members)
        coords[region] = (lat, lon)
        breakdown[region] = counts

    with_date = "--with-date" in sys.argv

    months = sorted({meta[n][0][:7] for n in names})
    trait_labels = list(LOCATION_ORDER) + months + (["Date"] if with_date else [])

    def row_values(name):
        date, loc = meta[name]
        vals = ["1" if loc == t else "0" for t in LOCATION_ORDER]
        vals += ["1" if date[:7] == m else "0" for m in months]
        if with_date:
            vals.append(date)
        return vals

    width = max(len(n) for n in names) + 3

    # --- NEXUS ---------------------------------------------------------------
    out = os.path.join(HERE, "Diamond_Princess.nex")
    with open(out, "w") as w:
        w.write("#nexus\n\n")
        w.write("[!\n")
        w.write("SARS-CoV-2 haplotypes from the Diamond Princess cruise-ship outbreak and\n")
        w.write("contemporaneous sequences, as distributed with:\n\n")
        w.write("Bo Xu, Lun Li, Cuiping Li, Anke Wang, Zhuojing Fan, Shuhui Song, HaploThread:\n")
        w.write("a scalable desktop tool for efficient haplotype network inference and\n")
        w.write("interactive visualization, Molecular Biology and Evolution, Volume 43, Issue 3,\n")
        w.write("March 2026, msag052, https://doi.org/10.1093/molbev/msag052\n\n")
        w.write("Sequences from Diamond_Princess.phy (%d taxa, %d variable sites).\n" % (ntax, nchar))
        w.write("Traits from Diamond_Princess.tsv; world coordinates added -- see README.md.\n")
        if with_date:
            w.write("Includes the categorical 'Date' trait: needs SplitsTree built after 2026-08-20.\n")
        else:
            w.write("The date is carried by the three numerical month-bin traits; run make-traits.py\n")
            w.write("--with-date for a categorical 'Date' trait with the exact value per sample.\n")
        w.write("]\n\n")

        w.write("BEGIN TAXA;\n")
        w.write("DIMENSIONS ntax=%d;\n" % ntax)
        w.write("TAXLABELS\n")
        for i, name in enumerate(names, start=1):
            w.write("\t[%d] '%s'\n" % (i, name))
        w.write(";\n")
        w.write("END; [TAXA]\n\n")

        w.write("BEGIN TRAITS;\n")
        w.write("DIMENSIONS [nTax=%d] nTraits=%d;\n" % (ntax, len(trait_labels)))
        w.write("FORMAT labels=yes missing=? separator=Comma;\n")
        w.write("[Coordinates apply per trait. The month bins%s carry 0/0 and are therefore\n"
                % (" and 'Date'" if with_date else ""))
        w.write(" skipped by the World Map view; the 8 locations come first so that all of\n")
        w.write(" them are drawn.]\n")
        w.write("TRAITLATITUDE")
        for t in trait_labels:
            w.write(" %s" % (("%.4f" % coords[t][0]) if t in coords else "0.0"))
        w.write(";\n")
        w.write("TRAITLONGITUDE")
        for t in trait_labels:
            w.write(" %s" % (("%.4f" % coords[t][1]) if t in coords else "0.0"))
        w.write(";\n")
        w.write("TRAITLABELS\n")
        w.write(" " + " ".join("'%s'" % t for t in trait_labels) + ";\n")
        w.write("MATRIX\n")
        for i, name in enumerate(names, start=1):
            w.write("[%d] '%s'%s%s\n" % (i, name, " " * (width - len(name)),
                                         ",".join(row_values(name))))
        w.write(";\n")
        w.write("END; [TRAITS]\n\n")

        w.write("BEGIN CHARACTERS;\n")
        w.write("DIMENSIONS ntax=%d nchar=%d;\n" % (ntax, nchar))
        w.write("FORMAT\n")
        w.write("\tdatatype='DNA' missing=? gap=- symbols=\"acgt\" labels=left transpose=no interleave=no;\n")
        w.write("MATRIX\n")
        for name in names:
            w.write("'%s'%s%s\n" % (name, " " * (width - len(name)), seqs[name].lower()))
        w.write(";\n")
        w.write("END; [CHARACTERS]\n")

    # --- the trait matrix ----------------------------------------------------
    #
    # A header line of trait names, a Coordinates line giving each trait its own
    # 'lat,long' pair, then one line of values per taxon. File -> Import -> Taxon Traits
    # has always read this; the Open Haplotype Data dialog does too, from 2026-08-20 --
    # both go through ImportTaxonTraits.parse. Per-trait coordinates are what that
    # dialog's older one-line-per-taxon layout could not express: there a line's single
    # (lat,long) applied to every trait named on it, so the month bins could not be left
    # unplaced, and they would have shown up on the world map.
    out2 = os.path.join(HERE, "Diamond_Princess-traits-matrix.tsv")
    with open(out2, "w") as w:
        w.write("traits\t" + "\t".join(trait_labels) + "\n")
        w.write("Coordinates\t" + "\t".join(
            ("%.4f,%.4f" % coords[t]) if t in coords else "0.0,0.0"
            for t in trait_labels) + "\n")
        for name in names:
            w.write(name + "\t" + "\t".join(row_values(name)) + "\n")

    # --- NEXUS for PopART ------------------------------------------------------
    #
    # PopART's parser differs from SplitsTree's: ntax may not appear in the CHARACTERS
    # block ("ntax only allowed in characters block along with new taxa"), and its TRAITS
    # block is its own dialect (Dimensions NTRAITS, TraitLabels, a Matrix of comma-separated
    # counts), written here exactly in the shape of PopART's example file: no comments,
    # nothing quoted, one space between a name and its data, and labels of letters, digits
    # and '_' only ('/', '-' and '.' of the strain names become '_'). The trait coordinates
    # are left out, since they only feed PopART's map view. PopART draws every trait of the
    # block into one pie, so the file carries one kind of trait, the eight locations; the
    # month bins are left out.
    def popart_label(name):
        return re.sub(r"[^A-Za-z0-9_]", "_", name)

    popart_names = [popart_label(n) for n in names]
    if len(set(popart_names)) != len(popart_names):
        sys.exit("replacing '/' by '_' makes two taxon labels coincide")
    out3 = os.path.join(HERE, "Diamond_Princess-popart.nex")
    with open(out3, "w") as w:
        w.write("#NEXUS\n\n")
        w.write("BEGIN TAXA;\n")
        w.write("DIMENSIONS NTAX=%d;\n" % ntax)
        w.write("TAXLABELS\n")
        for name in popart_names:
            w.write("%s\n" % name)
        w.write(";\n")
        w.write("END;\n\n")
        w.write("BEGIN CHARACTERS;\n")
        w.write("DIMENSIONS NCHAR=%d;\n" % nchar)
        w.write("FORMAT DATATYPE=DNA MISSING=? GAP=-;\n")
        w.write("MATRIX\n")
        for name, pname in zip(names, popart_names):
            w.write("%s %s\n" % (pname, seqs[name].upper()))
        w.write(";\n")
        w.write("END;\n\n")
        w.write("BEGIN TRAITS;\n")
        w.write("Dimensions NTRAITS=%d;\n" % len(LOCATION_ORDER))
        w.write("Format labels=yes missing=? separator=Comma;\n")
        w.write("TraitLabels " + " ".join(LOCATION_ORDER) + ";\n")
        w.write("Matrix\n")
        for name, pname in zip(names, popart_names):
            w.write("%s %s\n" % (pname, ",".join(row_values(name)[:len(LOCATION_ORDER)])))
        w.write(";\n")
        w.write("END;\n")

    # --- report --------------------------------------------------------------
    print("taxa %d, characters %d, traits %d" % (ntax, nchar, len(trait_labels)))
    print("\n%-16s %6s  %10s %11s" % ("location", "n", "latitude", "longitude"))
    for t in LOCATION_ORDER:
        n = sum(1 for x in names if meta[x][1] == t)
        print("%-16s %6d  %10.4f %11.4f" % (t, n, coords[t][0], coords[t][1]))
    print("\nmonth bins: " + ", ".join("%s (%d)" % (m, sum(1 for x in names if meta[x][0][:7] == m))
                                       for m in months))
    print("\nregional centroids are weighted means of:")
    for region in REGIONS:
        print("  %-14s %s" % (region, ", ".join("%s x%d" % (c, n)
                                                for c, n in breakdown[region].most_common())))
    print("\ncategorical 'Date' trait: %s" %
          ("included (--with-date)" if with_date else "omitted; pass --with-date to include it"))
    print("\nwrote:\n  %s\n  %s\n  %s" % (out, out2, out3))


if __name__ == "__main__":
    main()
