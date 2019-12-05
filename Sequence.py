#!/usr/bin/python3.8

import os
import random
from string import ascii_uppercase, digits

from Bio import Seq, SeqUtils, SeqIO, SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIWWW

__author__ = 'Milosz Chodkowski'
__license__ = "MIT"
__version__ = "1.0"
__status__ = "Production"


class Sequence:
    dna_bases = IUPAC.IUPACUnambiguousDNA.letters
    rna_bases = IUPAC.IUPACUnambiguousRNA.letters
    amino_acids = IUPAC.IUPACProtein.letters

    def __init__(self, size: int = 100, seq_type: str = 'D', id: str = None, seq=None):
        """Creates random Sequence of given size and type.

        :param size: Size of sequence.
        :param seq_type: Sequence type, D = DNA, R = RNA, P = Protein.
        :param id: ID of sequence.
        :param seq: Ready Sequence object.
        """
        self.s_type = {'D': 'DNA', 'R': 'RNA', 'P': 'PROTEIN'}[str(seq_type)]

        # If size is not multiple of 3, then make it bigger
        self.size = size if not size % 3 else size + (3 - size % 3)
        # If sequence is not none and if it's instance of Sequence class
        self.seq = seq if seq and isinstance(seq, Sequence) else self.generate_sequence()

        self.id = id if id else ''.join(random.choice(ascii_uppercase) for i in range(2)) + '|' \
                                + ''.join(random.choice(digits) for i in range(random.randint(4, 7)))

        self.record = SeqRecord.SeqRecord(self.seq, id=self.id)

    def show(self):
        print('Sequence: {}\nID: {}'.format(self.seq, self.id))

    def generate_sequence(self):
        """Generates random sequence based on type.

        :return: Bio.Seq object.
        """
        if self.s_type not in {'DNA', 'RNA', 'PROTEIN'}:
            raise ValueError('Wrong type of sequence')
        else:
            if self.s_type == 'DNA':
                seq = Seq.Seq(''.join(random.choice(Sequence.dna_bases) for i in range(self.size)))
            elif self.s_type == 'RNA':
                seq = Seq.Seq(''.join(random.choice(Sequence.rna_bases) for i in range(self.size)))
            else:
                seq = Seq.Seq(''.join(random.choice(Sequence.amino_acids) for i in range(self.size)))
        return seq

    def calculate_gc(self):
        """Calculates the GC percent in sequence.

        :return: Float number - GC percent.
        """
        if self.s_type == 'PROTEIN':
            raise TypeError('GC are not in {} sequence'.format(self.s_type))
        return SeqUtils.GC(self.seq)

    def transcribe(self):
        """Transcribes to RNA sequence if sequence is type D (DNA).

        :return: Seq object of type RNA.
        """
        if self.s_type != 'DNA':
            raise TypeError('Sequence type {} can not be transcribed.'.format(self.s_type))
        return Seq.Seq.transcribe(self.seq)

    def translate(self):
        """Translates to Protein sequence if sequence type is R (RNA).

        :return: Seq object of type Protein.
        """
        if self.s_type != 'RNA':
            raise TypeError('Sequence type {} can not be translated.'.format(self.s_type))
        return Seq.Seq.translate(self.seq)

    def reversed_transcription(self):
        """Given the seq of type RNA transcribes it to DNA.

        :return: Seq object of type DNA.
        """
        if self.s_type != 'RNA':
            raise TypeError('Sequence type {} can not be transcribed in reverse.'.format(self.s_type))
        return Seq.back_transcribe(self.seq)

    def get_sequence_elems(self):
        """Creates iterator of all bases in sequence.

        :return: Sequence bases iterator.
        """
        for base in self.seq:
            yield base

    def get_complement(self):
        """Gives the complement strand of sequence.

        :return: Complement Seq object.
        """
        return Seq.reverse_complement(self.seq)

    def save_to_fasta(self, fn=None, description='None'):
        """Saves sequence to file in fasta format.

        :param fn: Filename.
        :param description: Record description
        :return: None
        """
        if fn is None:
            fn = '{}.fasta'.format(self.record.id)
        self.record.description = description
        try:
            SeqIO.write(self.record, handle=fn, format='fasta')
        except OSError as exc:
            raise exc

    def read_from_fasta(self, fn=None):
        """Reads SeqRecord from file.

        If given file doesn't exists, the method takes first file in current directory.
        :param fn: Filename of fasta file.
        :return: True if file was loaded, else False
        """
        if fn is None:
            for file in os.listdir(os.curdir):
                if not file.endswith('.fasta'):
                    continue
                self.record = SeqIO.read(file, 'fasta')
                self.seq = self.record.seq
                self.id = self.record.id
                return True
        else:
            self.record = SeqIO.read(fn, 'fasta')
            self.seq = self.record.seq
            self.id = self.record.id
            return True
        return False

    def use_in_blast(self):
        if self.s_type == 'DNA' or self.s_type == 'RNA':
            income = NCBIWWW.qblast('blastn', 'nt', self.record.format('fasta'))
        elif self.s_type == 'PROTEIN':
            income = NCBIWWW.qblast('blastp', 'pdb', self.record.format('fasta'))
        else:
            raise TypeError('Type {} is not valid to use in blast'.format(self.s_type))
        result = income.read()
        return result
