tree = None

def build(node, node_l, node_r, arr):
    if node_l == node_r - 1:
        tree[node] = arr[node_l]
    else:
        mid = node_l + (node_r - node_l) // 2
        build(2 * node + 1, node_l, mid, arr)
        build(2 * node + 2, mid, node_r, arr)
        tree[node] = min(tree[2 * node + 1], tree[2 * node + 2])

def get(node, node_l, node_r, l, r): # [node_l; node_r), [l, r)
    if l <= node_l and node_r <= r: # вершина целиком внутри запроса
        return tree[node]

    if node_r <= l or r <= node_l: # вершина целиком снаружи запроса
        return int(1e18) # нейтральный элемент

    mid = node_l + (node_r - node_l) // 2
    # вершина 2 * node + 1 отвечает за [node_l, mid)
    # вершина 2 * node + 2 отвечает за [mod, node_r)

    return min(get(2 * node + 1, node_l, mid, l, r),
               get(2 * node + 2, mid, node_r, l, r))
    
def update(node, node_l, node_r, i, new_value):
    if node_l == node_r - 1: # лист
        assert i == node_l
        tree[node] = new_value
        return

    mid = node_l + (node_r - node_l) // 2

    if i < mid:
        update(2 * node + 1, node_l, mid, i, new_value)
    else:
        update(2 * node + 2, mid, node_r, i, new_value)
    tree[node] = min(tree[2 * node + 1], tree[2 * node + 2])
    
n = int(input())
tree = [0 for i in range(4*n)]
build(0, 0, n, [int(x) for x in input().split()])

while True:
    try:
        op, l, r = input().split()
        l = int(l) - 1
        r = int(r)
    except:
        break

    if op == 'min':
        print(get(0, 0, n, l, r))
    else:
        update(0, 0, n, l, r)
