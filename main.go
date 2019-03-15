package main

import (
	"bufio"
	"fmt"
	"os"
	"strings"
	"sync"
	"time"
	"github.com/andrewwhwang/go-radix"
	lev "github.com/andrewwhwang/levenshtein"
	"github.com/Workiva/go-datastructures/bitarray"
)

type pair struct {
    srPos, pos int
}

func min(int1, int2 int) int {
	if int1 < int2 {
		return int1
	}
	return int2
}

func max(int1, int2 int) int {
	if int1 > int2 {
		return int1
	}
	return int2
}

func reverse(s string) string {
    runes := []rune(s)
    for i, j := 0, len(runes)-1; i < j; i, j = i+1, j-1 {
        runes[i], runes[j] = runes[j], runes[i]
    }
    return string(runes)
}

//reads shortreads line by line
func getReadsFQ(filename string) chan string{

	file, _ := os.Open(filename)
	ch := make(chan string)
	scanner := bufio.NewScanner(file)

	go func() {
		defer close(ch)
		defer file.Close()
		for i := 0; scanner.Scan(); i++ {
			if i % 4 == 1 {
				ch <- scanner.Text()
			}
		}
	}()
	return ch
}

//make radix tree of short reads
func makeTree(filename string) *radix.Tree {
	tree := radix.New()
	for word := range getReadsFQ(filename) {
		counter := 1
		if v, i := tree.Get(word); i == true {
			counter = v.(int) + 1
		}
		tree.Insert(word, counter)
	}
	return tree
}

//Walks though leafs, outputting sorted shortreads
func uniqueIter(tree *radix.Tree) chan string {
	uniqueSRs := make(chan string)
	go func() {
		defer close(uniqueSRs)
		tree.Walk(func(s string, v interface{}) bool {
			uniqueSRs <- s
			// fmt.Println(v)
			// fmt.Print(s + " ")
			return false
		})
	}()
	return uniqueSRs
}

//make hashtable and bloom filter of kmer
func makeTables(filename string, k int) ([][]int, bitarray.BitArray, string) {
	
	arraySize := 1 << (2 * uint(k))
	bloom := bitarray.NewBitArray(uint64(arraySize))
	hashtable := make([][]int, arraySize, arraySize)

	file, _ := os.Open(filename)
	defer file.Close()

	scanner := bufio.NewScanner(file)
	scanner.Scan()
	scanner.Scan()
	ref := scanner.Text()

	counter := 0
	for hash := range rollingHash(&ref, k, 1) {
		hashtable[hash] = append(hashtable[hash], counter)
		bloom.SetBit(uint64(hash))
		counter++
	}

	return hashtable, bloom, ref
}

func baseMap(s byte) int {
	switch s {
	case 'A':
		return 0
	case 'C':
		return 1
	case 'G':
		return 2
	case 'T':
		return 3
	default:
		panic("Non-nucleotide Detected")
	}
}

func rollingHash(s *string, k, kmerInterval int) chan int {
	hashChan := make(chan int)
	arraySize := 1 << (2 * uint(k))
	rollHash := 0
	
	//get initial hash for substring(0, k-1]
	for i := 0; i < k ; i++ {
		rollHash += baseMap((*s)[i])
		rollHash <<= 2
	}
	
	//rolling window
	go func() {
		defer close(hashChan)
		for i := k; i < len(*s); i++ {
			rollHash += baseMap((*s)[i])
			rollHash &= arraySize - 1
			if (i % kmerInterval == 0) {
				hashChan <- rollHash
			}
			rollHash <<= 2
		}
	}()
	return hashChan
}
func extendTail(sr string, ref *string, srPos, refPos, k, thres, window int) int{
	end := srPos + k
	offset := refPos - srPos

	for i := srPos + k - window; i < len(sr) - window; i++{
		eDist := lev.Distance(sr[i:i+window], (*ref)[i+offset:i+offset+window])
		end++
		if eDist > thres || i + offset + window >= len(*ref) - 1 {
			break
		}
	}
	return end
}

func extendHead(sr string, ref *string, srPos, refPos, k, thres, window int) int{
	start := srPos
	offset := refPos - srPos

	for i := srPos; i >= 0; i--{
		eDist := lev.Distance(sr[i:i+window], (*ref)[i+offset:i+offset+window])
		start = i
		if eDist > thres || i + offset <= 0 {
			break
		}
	}
	return start
}

// func extendHead(sr string, ref *string, srPos, refPos, k, thres, window int) int{
// 	start := srPos
// 	eList := make([]int, window)
// 	var eChan chan int

// 	wordLen := min(srPos, refPos)
// 	srCrop := reverse(sr[srPos-wordLen:srPos])
// 	refCrop := reverse((*ref)[refPos-wordLen:refPos])

// 	if wordLen < 32 {
// 		eChan = lev.MyerDist(srCrop, refCrop)
// 	} else {
// 		eChan = lev.MyerDistDiag(srCrop, refCrop, 15)
// 	}

// 	for wordLen > 0 && eList[len(eList)-1] - eList[0] <= thres{
// 		start--
// 		wordLen--
// 		eDist := <-eChan
// 		eList = append(eList[1:], eDist)
// 	}
// 	return start
// }

// func extendTail(sr string, ref *string, srPos, refPos, k, thres, window int) int{
// 	end := srPos + k
// 	eList := make([]int, window)
// 	var eChan chan int

// 	wordLen := min(len(sr) - (srPos + k), len(*ref) - (refPos + k))
// 	srCrop := sr[srPos+k:srPos+k+wordLen]
// 	refCrop := (*ref)[refPos+k : refPos+k+wordLen]

// 	if wordLen < 64 {
// 		eChan = lev.MyerDist(srCrop, refCrop)
// 	} else {
// 		eChan = lev.MyerDistDiag(srCrop, refCrop, 15)
// 	}

// 	for wordLen > 0 && eList[len(eList)-1] - eList[0] <= thres{
// 		end++
// 		wordLen--
// 		eDist := <-eChan
// 		eList = append(eList[1:], eDist)
// 	}
// 	// close(eChan)
// 	return end
// }

//TODO: make the diagonal width offset dynamic
// make a square matrix from ranging from [0 - sr head] and [sr tail - end]
// rolling average window along the diagonal, break when the distance spikes
func getFuzzy(sr string, ref *string, srPos , refPos, k, window, thres int) (string, int) {
	
	if window > k {
		panic("window larger than k")
	}

	offset := refPos - srPos
	// refLen := len(*ref)

	//fuzzy extension going left
	var wg sync.WaitGroup
	var start, end int
	wg.Add(2)

	//fuzzy extension going left
	go func() {
		start = extendHead(sr, ref, srPos, refPos, k, thres, window)
		wg.Done()
	}()
	
	//fuzzy extension going right
	go func() {
		end = extendTail(sr, ref, srPos, refPos, k, thres, window)
		wg.Done()
	}()

	wg.Wait()
	return sr[start:end], start+offset
}


func printAln(aln string, offset int) {
	if offset > 0 {	
		padding := strings.Repeat(" ", offset)
		fmt.Println(padding + aln)
	} else {
		fmt.Println(aln[-1 * offset:])
	}
}


func main() {
	start := time.Now()
	k, kmerInterval := 8, 4
	window, thres := 6, 3

	//make radix tree of short reads
	tree := makeTree("resources/short_reads.fq")
		
	//make k-mer dictionary of reference reads. buckethash + bloom filter
	kmerHash, bloom, ref := makeTables("resources/ref.fa", k)

	fmt.Println(ref)
	// //interate through unique short reads and align them

	for sr := range uniqueIter(tree) {
		go func(sr string) {
			var srPos, bestAlnLen, bestOffset int
			var bestAln string

			for hash := range rollingHash(&sr, k, kmerInterval) {
				if in, _ := bloom.GetBit(uint64(hash)); !in{ 
					srPos += kmerInterval
					continue 
				}

				positions := kmerHash[hash]
				// fmt.Println(positions)
				//skip if too many kmers
				if len(positions) < 10 {
					for _, pos := range positions {
						aln, offset := getFuzzy(sr, &ref, srPos, pos, k, window, thres)
						alnLen := len(aln)
						if alnLen > bestAlnLen {
							bestAln, bestAlnLen, bestOffset = aln, alnLen, offset
						}
					}
				}
				srPos += kmerInterval
			}
			if bestAlnLen > 10 {
				//do smith-waterman alignment
				printAln(bestAln, bestOffset)
			}
		}(sr)
	}    
	elapsed := time.Since(start)
    fmt.Println(elapsed)
	fmt.Println("done")
}