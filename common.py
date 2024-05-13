#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2024 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

from ssw import AlignmentMgr

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
