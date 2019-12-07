#!/usr/bin/env python3
"""
Plot module contains methods for dot-plotting sequences.
"""
import sys

import matplotlib.pylab as pylab
from Bio import SeqIO
from Bio.Blast import NCBIXML


def plot_from_xml(fn1=None, fnx=None, window=7):
    """Makes a dot-plot from xml and fasta sequence.

    :param fn1: Fasta file.
    :param fnx: XML file.
    :param window: Threshold, for example (5)'
    :return: None
    """
    if fn1 is not None and fnx is not None:
        with open(fn1, 'r') as fh, open(fnx, 'r') as fx:
            window = int(window)
            file1 = NCBIXML.read(fh)
            filex = NCBIXML.read(fx)
            seq_one = str(file1.seq).upper()
            seq_two = str(filex.seq).upper()
            data = [[(seq_one[i:i + window] != seq_two[j:j + 5]) for j in range(len(seq_one) - window)] for i
                    in range(len(seq_two) - window)]

            pylab.gray()
            pylab.imshow(data)
            pylab.xlabel(
                '{} (length {} bp)'.format(file1.seq, len(file1.seq)))
            pylab.ylabel(
                '{} (length {} bp)'.format(filex.seq, len(filex.seq)))
            pylab.title('Dot plot using window size {}\n(allowing no mis-matches)'.format(window))
            pylab.show()


def plot_from_fasta(fn1=None, fn2=None, window=7):
    """Makes a dot-plot two fasta files.

    :param fn1: Fasta file.
    :param fn2: Fasta file.
    :param window: Threshold, for example (5)
    :return: None
    """
    if fn1 is not None and fn2 is not None:
        with open(fn1, 'r') as fh1, open(fn2, 'r') as fh2:
            window = int(window)
            file1 = SeqIO.read(fh1, format='fasta')
            file2 = SeqIO.read(fh2, format='fasta')
            dict_one = {}
            dict_two = {}
            for (seq, section_dict) in [(str(file1.seq).upper(), dict_one),
                                        (str(file2.seq).upper(), dict_two)]:
                for i in range(len(seq) - window):
                    section = seq[i:i + window]
                    try:
                        section_dict[section].append(i)
                    except KeyError:
                        section_dict[section] = [i]
            matches = set(dict_one).intersection(dict_two)
            print('{} unique matches'.format(len(matches)))
            x, y = [], []
            for section in matches:
                for i in dict_one[section]:
                    for j in dict_two[section]:
                        x.append(i)
                        y.append(j)
            pylab.cla()
            pylab.gray()
            pylab.scatter(x, y)
            pylab.xlim(0, len(file1) - window)
            pylab.ylim(0, len(file2) - window)
            pylab.xlabel('{} (length {} bp)'.format(file1.id, len(file1.seq)))
            pylab.ylabel('{} (length {} bp)'.format(file2.id, len(file2.seq)))
            pylab.title('Dot plot using window size {}\n(allowing no mis-matches)'.format(window))
            pylab.show()


if __name__ == '__main__':
    if sys.argv[4] == 'xml':
        plot_from_xml(sys.argv[1], sys.argv[2], sys.argv[3])
    elif sys.argv[4] == 'fasta':
        plot_from_fasta(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        pass
