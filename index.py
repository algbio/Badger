#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2024 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

from collections import defaultdict
import math

class QGramIndex:
    
    RANK = {'A' : 0,
            'C' : 1,
            'G' : 2,
            'T' : 3}
    
    def __init__(self, threshold, bc_len, q = 2):
        self.q = q
        # print("k:", self.q)
        self.threshold = bc_len - q + 1 - q*threshold
        if self.threshold <= 0:
            self.threshold = 4
        self.index = []
        self.list_of_q_grams() 

    def add_to_index(self, barcode, number):
        kmer = self.rank(barcode[:self.q])
        self.index[kmer][number] += 1
        for i in range(self.q, len(barcode)):
            kmer = self.update_rank(kmer, barcode[i])
            self.index[kmer][number] += 1
        return
        
    def list_of_q_grams(self):
        for i in range(0, pow(4,self.q)):
            self.index.append(defaultdict(int))
        return
        
    def get_closest(self, barcode, number):
        distances = defaultdict(int)
        kmers = []
        kmer = self.rank(barcode[:self.q])
        kmers.append(kmer)
        for i in range(self.q+1, len(barcode)):
            kmer = self.update_rank(kmer, barcode[i])
            kmers.append(kmer)
        options = set()
        for i in kmers:
            options.update(self.index[i].keys())
        return options
        for i in options:
            for j in range(0, pow(4,self.q)):
                x = abs(self.index[j].get(number,0) - self.index[j].get(i,0))
                distances[i] += x
        
        results = []
        for i in range(0, len(self.index[0])):
            if distances[i] <= self.threshold:
                if i != number:
                    results.append(i)
        
        return results
    
    def rank(self, qgram):
        rank = 0
        for i in range(0,self.q):
            rank += QGramIndex.RANK[qgram[i]]*(pow(4,i))
        return rank
    
    def update_rank(self, rank, b):
        return math.floor(rank/4) + QGramIndex.RANK[b] * (pow(4,self.q-1))
        
    def get_close(self, barcode, number):
        distances = defaultdict(int)
        kmers = []
        kmer = self.rank(barcode[:self.q])
        kmers.append(kmer)
        for i in range(self.q, len(barcode)):
            kmer = self.update_rank(kmer, barcode[i])
            kmers.append(kmer)
        for i in kmers:
            for j in self.index[i].keys():
                if j > number:
                    distances[j] += self.index[i][j]
        results = []
        for bc in distances.keys():
            if distances[bc] >= self.threshold:
                results.append(bc)
        return results

class KMerIndex:
    
    def __init__(self, k):
        self.k = k
        self.index = defaultdict(list)
        
    def _get_kmers(self, seq):
        if len(seq) < self.k:
            return
        kmer = seq[:self.k]
        yield kmer
        for i in range(self.k, len(seq)):
            kmer = kmer[1:] + seq[i]
            yield kmer
            
    def add_to_index(num, barcode):
        for kmer in self._get_kmers(barcode):
            self.index[kmer].append(num)
            
    def get_occurrences(self, sequence, max_hits=0, min_kmers=1, hits_delta=1, ignore_equal=False):
        barcode_counts = defaultdict(int)
        barcode_positions = defaultdict(list)
        for pos, kmer in enumerate(self._get_kmers(sequence)):
            for i in self.index[kmer]:
                barcode_counts[i] += 1
                barcode_positions[i].append(pos)

        result = []
        for i in barcode_counts.keys():
            count = barcode_counts[i]
            if count < min_kmers:
                continue
            if ignore_equal and self.seq_list[i] == sequence:
                continue
            result.append((self.seq_list[i], count, barcode_positions[i]))

        return result
