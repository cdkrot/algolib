def siftUp(cur):
	while cur > 1 and H[cur] < H[cur // 2]:
		swap(H[cur], H[cur // 2])
		cur //= 2

def siftDown(cur):
	while cur * 2 <= sz:
		mv = cur * 2
		if cur * 2 + 1 <= sz and H[cur * 2 + 1] < H[mv]:
			mv = cur * 2 + 1
		if H[cur] <= H[mv]:
			break
		swap(H[cur], H[mv])
		cur = mv

def getMin():
	return H[1]

def add(x):
	sz += 1
	H[sz] = x
	siftUp(sz)

def extractMin():
	swap(H[1], H[sz])
	sz -= 1
	siftDown(1)
	return H[sz + 1]
