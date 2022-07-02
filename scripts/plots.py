import os
import sys
import argparse
from pathlib import Path

import pandas
import matplotlib.pyplot as plt


def parse_seqtkfqchk(infile):

    with open(infile, 'r') as f:
        lines = f.read().splitlines()

    basic_stats = {}
    read_summary = {}
    histo = []
    for l in lines:
        if l.startswith("min_len"):
            pairs = l.split(';')
            pairs.pop()
            for x in pairs:
                stat, value = x.split(':')
                basic_stats[stat] = value
        elif l.startswith("POS"):
            keys = l.split()
        elif l.startswith("ALL"):
            values = l.split()
            for k, v in zip(keys, values):
                basic_stats[k] = v
        elif l[:1].isdigit():
            histo.append(l.split())

    return basic_stats, histo


def plot_lengths(seqtkinfile):
    flowcell = os.path.basename(Path(seqtkinfile)).rstrip("_seqtkfqchk.txt")
    basic_stats, histo = parse_seqtkfqchk(seqtkinfile)
    total = basic_stats['#bases']

    df = pandas.DataFrame(histo, columns=list(basic_stats.keys())[3:])
    df['#bases'] = df["#bases"].astype(int)
    df['#bases'].plot()

    plt.title("Readlength")
    plt.legend([flowcell], loc=1, fontsize='x-small', bbox_to_anchor=(1, 1))
    outfile = f"{seqtkinfile}.png"
    plt.savefig(outfile)


def get_parser():
    """Return an argparse instance"""

    p = argparse.ArgumentParser("Parse and plot seqtq fqchk output.")
    p.add_argument('--debug', action='store_true',
                   help="Print debug logging to stdout")
    p.add_argument('seqtkfile', type=str, help="path to seqtk fqchk output.")

    return p


def main(args):
    parser = get_parser()
    args = parser.parse_args()
    plot_lengths(args.seqtkfile)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
