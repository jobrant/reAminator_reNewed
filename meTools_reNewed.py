#!/usr/bin/env python3
# (c) 2013, Russell Darst, University of Florida
# Optimized version

from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
from typing import List, Dict, Tuple, Union
from dataclasses import dataclass
from functools import lru_cache

@dataclass
class Patch:
    A: float
    B: float
    C: float
    D: float
    V: int

def quilt(pattern: str) -> List[Patch]:
    if not pattern.strip():
        return []

    P = pattern + ' '
    Q = [Patch(np.nan, i, i + P[i:].find(' ') - 1, np.nan, '#*'.find(j))
         for i, j in enumerate(P) if ' ' == P[i-1] != j]
    Q = [Patch(np.nan, np.nan, np.nan, np.nan, -1)] + Q + [Patch(np.nan, np.nan, np.nan, np.nan, -1)]

    for i in range(1, len(Q)):
        Q[i].A = Q[i-1].C
        Q[i-1].D = Q[i].B

    return Q[1:-1]

def read_FASTA(which_file: str, **kwargs) -> 'meTable':
    with open(which_file, 'r') as handle:
        seqs = list(SeqIO.parse(handle, 'fasta'))
    return meTable(seqs, **kwargs)

def seek(obj: Union[str, 'meTable'], col: Union[int, str] = -1, keep: str = 'True', drop: str = 'False'):
    keep_func = eval(f'lambda v, ad, bc: {keep}')
    drop_func = eval(f'lambda v, ad, bc: {drop}')
    
    if isinstance(obj, str): 
        obj = read_FASTA(obj)
    if isinstance(col, str): 
        col = obj.__head__.index(col)

    test = lambda f, x: f(x.V, x.D - x.A - 1, x.C - x.B + 1)
    for seq in obj:
        patt = quilt(seq[col])
        if not patt:
            continue

        q = patt[0]
        for p in patt[1:]:
            if test(drop_func, p) or p.V == q.V:
                q.C, q.D = p.C, p.D
                continue

            if test(keep_func, q):
                ends = slice(q.B, q.C + 1)
                yield (seq[0], seq[1], q.A, q.B, q.C, q.D,
                       seq[col][ends].count('#'), seq[col][ends].count('*'))

            q = p

class meTable:
    def __init__(self, contig: List[SeqIO.SeqRecord], BS: float = 95, unambig: bool = True, **kwargs):
        self.BS = BS / 100
        self.__head__ = ['locus', 'seqID']
        self.__offsets__ = {}
        self.__sites__ = {}
        self.ref = contig[0]
        self.seqs = {seq.id: seq for seq in contig[1:]}

        self.__rules__ = {
            'locus': lambda x: self.ref.id,
            'seqID': lambda x: x.id
        }

        if not kwargs:
            kwargs = dict(CG=1, GC=2)
        for k, v in kwargs.items():
            if not isinstance(v, int) or v == 0:
                raise ValueError('Site parameters must be in format GC=2, etc.')

            self.__head__.append(k)
            self.__sites__[k] = self.__match__(k, v)
            self.__offsets__[k] = v
            self.__rules__[k] = lambda x, y=k: self.__score__(x, y)

        self.__head__.sort()
        self.__checkBS__()

        if unambig:
            self.__sites__ = {k: self.__pare__(k, self.__sites__[k])
                              for k in kwargs.keys()}

    def __checkBS__(self):
        B = self.BS / (1 - self.BS)
        C = self.__pare__('C', self.__match__('C', 1))
        D = lambda x, y: sum(x[c].upper() == y for c in C)
        self.seqs = {i: j for i, j in self.seqs.items() if D(j, 'T') >= B * D(j, 'C')}

    def __getitem__(self, i):
        if isinstance(i, str):
            return [self.__rules__[j](self.seqs[i]) for j in self.__head__]
        else:
            return [self.__rules__[j](list(self.seqs.values())[i])
                    for j in self.__head__]

    def __len__(self):
        return len(self.seqs)

    @lru_cache(maxsize=None)
    def __match__(self, M: str, N: int) -> List[int]:
        ref = self.ref.seq.upper()
        return [n + N - 1 for n in range(len(ref)) if ref.find(M, n) == n]

    def __pare__(self, K: str, I: List[int]) -> List[int]:
        return sorted(set(I) - set(
            j for k in self.__head__[2:] if k != K for j in self.__sites__[k]))

    def __score__(self, read: SeqIO.SeqRecord, site: str) -> str:
        check = lambda x: (x == 'C') - (x == 'T')
        ref_len, seq = len(self.ref.seq), read.seq
        I = self.__sites__[site]
        J = [check(seq[i].upper()) for i in I]
        K = [i for i, j in enumerate(J) if j != 0]
        if not K:
            return ' ' * ref_len

        L = [' ' * I[K[0]]]
        for i, j in zip(K[:-1], K[1:]):
            L.append(' *#'[J[i]] + ' +-'[(J[i] + J[j]) // 2] * (I[j] - I[i] - 1))
        L.append(' *#'[J[-1]] + ' ' * (ref_len - I[K[-1]] - 1))
        return ''.join(L)

    def screen(self, i: str):
        n, seq = 0, self.seqs[i]
        while n < len(self.ref):            
            print('ref     ' + self.ref.seq[n:n+70])
            print('read    ' + seq.seq[n:n+70])
            for site in self.__head__[2:]:
                print(f'{site:<8}' + self.__rules__[site](seq)[n:n+70])
            n += 78
            print()

    def write(self, OUT: Union[str, 'TextIO'], header: bool = False):
        if isinstance(OUT, str): 
            OUT = open(OUT, 'w')
        OUT.write('\n'.join(['\t'.join(self.__head__)] * header + [
            '\t'.join(map(str, j)) for j in self[:]]))