#python implementation
from pygtrie import CharTrie
from bitarray import bitarray
from skbio.alignment import StripedSmithWaterman
import numpy as np
import editdistance

BASE_MAP = {'A':0,'C':1,'G':2,'T':3}

def fastaSeq(filename):
    with open(filename, "r") as f:
        next(f)
        return f.readline()
        
def makeTree(filename):
    t = CharTrie()
    t.enable_sorting()
    with open(filename) as f:
        next(f)
        for i, f in enumerate(f.readlines()):
            if f and i % 4 == 0:
                t[f.strip()] = t.setdefault(f.strip(), 0) + 1
    return t

class RollingHash:
    def __init__(self, string):
        self.k = len(string)
        h = 0
        for char in string:
            h <<= 2
            h += BASE_MAP[char]
        self.curhash = h

    # Returns the current hash value.
    def hash(self):
        return self.curhash

    # Updates the hash by removing previtm and adding nextitm. Returns the updated
    # hash value.
    def slide(self, nextitm):
        
        self.curhash <<= 2
        self.curhash += BASE_MAP[nextitm]
        self.curhash &= (4**self.k) - 1

        return self.curhash

# Create hash table of all kmer occurances, and bloom filter
def makeHashtable(ref, k):
    length = 4**k
    bloom = bitarray(length)
    bloom.setall(False)
    hashtable = np.empty((length,), dtype=object)
    hashtable[:] = [[] for _ in range(length)]

    refHash = RollingHash(ref[:k])
    hashtable[refHash.hash()].append(0)

    for idx in range(k,len(ref)):
        currHash = refHash.hash()
        refHash.slide(ref[idx])
        hashtable[currHash].append(idx-k)
        bloom[currHash] = True

    return hashtable, bloom

def extendHead(sr, ref, srPos, refPos, k, window, thres):
    start = srPos
    offset = refPos - srPos

    for i in range(srPos, 0, -1):
        eDist = editdistance.eval(sr[i:i+window], ref[i+offset:i+offset+window])
        start = i
        if eDist > thres or i + offset <= 0:
            break

    return start

def extendTail(sr, ref, srPos, refPos, k, window, thres):
    end = srPos + k
    offset = refPos - srPos

    for i in range(srPos+k-window, len(sr)-window):
        eDist = editdistance.eval(sr[i:i+window], ref[i+offset:i+offset+window])
        end += 1
        if eDist > thres or i + offset + window >= len(ref) - 1:
            break

    return end

def getFuzzy(sr, ref, srPos, refPos, k, window, thres):
    offset = refPos - srPos

    start, end = 0, 0

	#fuzzy extension going left and right
    start = extendHead(sr, ref, srPos, refPos, k, window, thres)
    end = extendTail(sr, ref, srPos, refPos, k, window, thres)

    return sr[start:end], start+offset

def printAln(aln, offset):
    if offset > 0:
        padding = " " * offset
        print(padding + aln)
    else:
        print(aln[-1 * offset:])

def hasSimilarNeighbor(tree, str1):
    # arbitrarily set threshold to 95% to be considered similar
    simThres = int(0.95 * len(str1))
    prefix = str1[:simThres]
    if len(list(tree.items(prefix=prefix))) > 1:
        return True
    return False

def smithWaterson(bestAln, ref):
    aln,offset,_,_ = bestAln
    targetStart = max(0, offset)
    target = ref[targetStart:targetStart + len(aln)]

    query = StripedSmithWaterman(target)
    alignment = query(aln)
    return alignment.aligned_query_sequence

if __name__ == '__main__':
    refFile = "resources/ref.fa"
    srFile = "resources/short_reads.fq"
    k = 8
    window = 5
    thres = 3

    assert k >= window and window >= thres

    ref = fastaSeq(refFile)
    print(ref)

    tree = makeTree(srFile)
    hashtable, bloom = makeHashtable(ref, k)
    prevAlns = []
    prevSimilar = False

    for sr in tree:
        if prevSimilar and hasSimilarNeighbor(tree, sr):
            alns = []
            for _, _,srPos,pos in prevAlns:
                aln, offset = getFuzzy(sr, ref, srPos, pos, k, window, thres)
                if len(aln) > len(sr) * 0.10:
                    alns.append((aln,offset,srPos,pos))
            if alns:
                bestAln = max(alns, key=lambda x:len(x[0]))
                alignedSR = smithWaterson(bestAln, ref)
                printAln(alignedSR, bestAln[1])
        else:
            prevAlns = []
            srHash = RollingHash(sr[:k])
            for srPos in range(len(sr) - k):
                if srPos % 4 == 0:
                    currHash = srHash.hash()
                    if bloom[currHash] == False or len(hashtable[currHash]) > 10:
                        srHash.slide(sr[srPos + k])
                        continue
                    for pos in hashtable[currHash]:
                        aln, offset = getFuzzy(sr, ref, srPos, pos, k, window, thres)
                        if len(aln) > len(sr) * 0.10:
                            prevAlns.append((aln,offset,srPos,pos))
                srHash.slide(sr[srPos + k])
            if prevAlns:
                bestAln = max(prevAlns, key=lambda x:len(x[0]))
                alignedSR = smithWaterson(bestAln, ref)
                printAln(alignedSR, bestAln[1])
        prevSimilar = hasSimilarNeighbor(tree, sr)


# ref.fa(700bp) + short_reads.fq(30,391 SRs @ ~450bp/SR) = 2min 11secs
