#!/usr/bin/env python3

import sys
import Plot
from Sequence import Sequence

__author__ = 'Milosz Chodkowski'
__license__ = "MIT"
__version__ = "1.0"
__status__ = "Production"


def main():
    sq = Sequence(size=5000, seq_type='D')
    sq.read_from_fasta('dr.fasta')
    sq.blast_search()
    Plot.plot_from_fasta('dr.fasta', 'dy.fasta', sys.argv[1])


if __name__ == '__main__':
    main()
