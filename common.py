#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2024 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

from ssw import AlignmentMgr

RANK = {'A' : 0,
        'C' : 1,
        'G' : 2,
        'T' : 3}

UNRANK = {0 : 'A',
		  1 : 'C',
		  2 : 'G',
		  3 : 'T'}
        
def rank(seq, length):
    rank = 0
    for i in range(0,length):
        rank += RANK[seq[i]]*(pow(4,i))
    return rank
    
def unrank(rk, length):
	
	seq = ""
	
	for i in range(length):
		
		c = rk % 4 
		
		seq = seq + UNRANK[c]
		rk = rk // 4
	
	return seq
	

def dfs(visited, node, component, edges):
    if not visited[node]:
        component.append(node)
        visited[node] = True
        for neighbor in edges[node]:
            if not visited[neighbor]:
                dfs(visited, neighbor, component, edges)
    return component
    
def dfs_without_recursion(visited, node, component, edges):
    stack = []
    stack.append(node)
    while stack:
        node = stack.pop()
        if not visited[node]:
            component.append(node)
            visited[node] = True
            for i in edges[node]:
                stack.append(i)
    return component, visited

def get_score(bc1, bc2):
    aln = AlignmentMgr()
    align_mgr = AlignmentMgr(match_score=3, mismatch_penalty=3)
    align_mgr.set_read(bc1)
    align_mgr.set_reference(bc2)
    alignment = align_mgr.align(gap_open=2, gap_extension=2)
    return alignment.optimal_score
