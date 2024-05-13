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
        print("k:", self.q)
        # instead of distance calculate number of matching k-mers
        # then if we have a certain number of them give it back
        # take minimum of occurences and add to the sum
        # self.threshold = 2*self.q*threshold
        self.threshold = bc_len - q + 1 - q*threshold
        if self.threshold <= 0:
            self.threshold = 4
        self.index = []
        self.list_of_q_grams() 
        # dict mapping all q-grams/k-mers to the barcodes it is in (+ how many times! don't know how to do this yet) important: k-mers are ranked so displayed as numbers
        # Actually this should probably be a list if I rank my q-grams
        # two possibilities: List of dicts, barcode: number of occurences or List of lists, list containing one entry for each occurence
        # I think with get and default 0 dict might be faster but need to check
        # collections defaultdict! Not sure if better for my possibilities tho
        # Gut: Ich muss nicht initialisieren
        # Schlecht: Ich bin mir nicht sicher wie ich am besten da drüber iteriere oder ob das mit dem += 1 dann funktioniert. Wenn ja dann vielleicht tauschen. Andrey fragen.
        
        # What if I also number the barcodes? I can number them by when I add them to the nodes and then I can just access by list position
        # entry in the list is amount of times the q-gram appears in the barcode, default 0
        # Actually no, I don't know the number of barcodes when I init
        # Maybe if I append a 0 every time in add_to_index?
        # That seems to be the same as the initialisation step in q-gram distance but I need way more of them because I keep a profile for each barcode, not just one

    def add_to_index(self, barcode, number):
        # for i in range(0, pow(4,self.q)):
            # self.index[i].append(0)
        kmer = self.rank(barcode[:self.q])
        self.index[kmer][number] += 1
        for i in range(self.q, len(barcode)):
            kmer = self.update_rank(kmer, barcode[i])
            self.index[kmer][number] += 1
            # if barcode == "ATGCTTAAGCACCACA":
                # print(barcode[:self.q])
                # print(kmer, self.index[kmer][number])
        #print(self.index)
        return
        
    
    def list_of_q_grams(self):
        for i in range(0, pow(4,self.q)):
            #self.index.append([])
            self.index.append(defaultdict(int))
        #print(self.index)
        return
        
    def get_closest(self, barcode, number):
        distances = defaultdict(int)
        # switch loops and choose only ones that share at least one k-mer
        # maybe look up which k-mers are in the barcodes by ranking again
        # expanding on first line: what if i go over it twice? Once kinda what Andrey does and and then the ones I keep with that and do mine?
        kmers = []
        kmer = self.rank(barcode[:self.q])
        kmers.append(kmer)
        for i in range(self.q+1, len(barcode)):
            kmer = self.update_rank(kmer, barcode[i])
            kmers.append(kmer)
        # options probably the issue
        # maybe increase k
        # compare with and without q-gram distance to make sure the distance still changes anything
        # otherwise maybe find another threshold
        options = set()
        #print(self.index)
        for i in kmers:
            options.update(self.index[i].keys())
            # for j in range(0, len(self.index[i])):
                # if self.index[i][j] != 0:
                    # options.add(j)
        return options
        #print(len(options) == len(self.index[kmer].keys()))
        for i in options:
            for j in range(0, pow(4,self.q)):
                x = abs(self.index[j].get(number,0) - self.index[j].get(i,0))
                distances[i] += x
        # for i in range(0, pow(4,self.q)):
            # for j in range(0, len(self.index[i])):
                # x = abs(self.index[i][number] - self.index[i][j])
                # distances[j] += x
        
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
        # bitwise shift or integer division
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
                distances[j] += self.index[i][j]
        results = []
        for bc in distances.keys():
            if distances[bc] >= self.threshold:
                #print(distances[bc])
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
