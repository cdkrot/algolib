#!/usr/bin/python3
# Alice [cdkrot.me] Sayutina (2022)

import sys
from queue import SimpleQueue

class Edge:
    def __init__(self, to, flow, cap, nxt):
        self.to = to
        self.flow = flow
        self.cap = cap
        self.nxt = nxt

def main():
    n, m = map(int, input().split())
    S = 0
    T = n - 1

    edges = []
    head = [-1 for i in range(n)]

    dist = None
    
    for _ in range(m):
        a, b, cap = map(int, input().split())
        a -= 1
        b -= 1

        edges.append(Edge(b, 0, cap, head[a]))
        head[a] = len(edges) - 1

        edges.append(Edge(a, 0, 0, head[b]))
        head[b] = len(edges) - 1

    def bfs():
        q = SimpleQueue()
        q.put(S)

        nonlocal dist
        dist = [None for v in range(n)]
        dist[S] = 0

        while not q.empty():
            v = q.get()

            e = head[v]
            while e != -1:
                if edges[e].flow < edges[e].cap and dist[edges[e].to] is None:
                    dist[edges[e].to] = dist[v] + 1
                    q.put(edges[e].to)
                e = edges[e].nxt

        return dist[T] is not None
        
    def find_flow(v, max_cap, head):
        if v == T:
            return max_cap

        while True:
            e = head[v]
            if e == -1:
                return 0
            
            limit = min(edges[e].cap - edges[e].flow, max_cap)

            if limit and dist[edges[e].to] == dist[v] + 1 \
            and (result := find_flow(edges[e].to, limit, head)):
                edges[e].flow += result
                edges[e^1].flow -= result
                return result
            else:
                head[v] = edges[e].nxt

    ans = 0
    # bfs() computes dist[] and returns True, if dist[T] != inf.
    while bfs():
        headcopy = list(head) # копировать массив head.

        while (delta := find_flow(S, int(1e9), headcopy)):
            ans += delta

    print(ans)
    for e in range(0, len(edges), 2):
        print(edges[e].flow)

main()
