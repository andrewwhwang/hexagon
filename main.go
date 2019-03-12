package main

import (
	"bufio"
	"fmt"
	"os"
	"strings"
	"sync"
	"time"
	"github.com/andrewwhwang/go-radix"
	"github.com/andrewwhwang/levenshtein"
	"github.com/Workiva/go-datastructures/bitarray"
)

type pair struct {
    srPos, pos int
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

//string distance distance rolling window
//https://github.com/xrash/smetrics ukkonen's 
//TODO: make the diagonal width offset dynamic
// make a square matrix from ranging from [0 - sr head] and [sr tail - end]
// rolling average window along the diagonal, break when the distance spikes
func getFuzzy(sr, ref *string, srPos , refPos, k, window, thres int) (string, int) {
	var start, end int

	if window > k {
		panic("window larger than k")
	}

	offset := refPos - srPos
	refLen := len(*ref)

	//fuzzy extension going left
	var wg sync.WaitGroup
	wg.Add(2)
	go func() {
		defer wg.Done()
		for i := srPos + window; i >= window; i-- {
			start = i-window
			refStart := start+offset
			refEnd := i+offset
			if refStart > 0 && thres <= levenshtein.MyerDist((*sr)[start:i], (*ref)[refStart:refEnd]) {
				break
			}
		}
	}()
	
	//fuzzy extension going right
	go func() {
		defer wg.Done()
		for i := srPos + k - window; i <= len(*sr) - window; i++ {
			end = i+window
			refStart := i+offset
			refEnd := end+offset
			if refEnd < refLen && thres <= levenshtein.MyerDist((*sr)[i:end], (*ref)[refStart:refEnd]){
				break
			}
		}
	}()

	wg.Wait()
	return (*sr)[start:end], start+offset
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
	window, thres := 8, 5

	//make radix tree of short reads
	//TODO: parallelize
	tree := makeTree("resources/truncated.fq")
		
	//make k-mer dictionary of reference reads. buckethash + bloom filter
	//TODO: parallelize kmerhash by sorting into 4 regions of table
	kmerHash, bloom, ref := makeTables("resources/ref.fa", k)

	fmt.Println(ref)
	// //interate through unique short reads and align them
	var prevCandidates []pair
	var prevSimilar bool
	for sr := range uniqueIter(tree) {
		var srPos, bestAlnLen, bestOffset int
		var bestAln string
		//if the new SR < 5% different, reuse the previous SR's seeds.
		diff := float32(tree.SuffixDifference(sr)) / float32(len(sr))
		if diff < 0.05 && prevSimilar {
			for _, prev := range prevCandidates {
				//check if seed positions are not out of index
				if prev.srPos+k < len(sr){
					aln, offset := getFuzzy(&sr, &ref, prev.srPos, prev.pos, k, window, thres)
					alnLen := len(aln)
					if alnLen > bestAlnLen {
						bestAln, bestAlnLen, bestOffset = aln, alnLen, offset
					}
				}
			}
		} else {
			
			//only run through process if > 5% different from last SR
			for hash := range rollingHash(&sr, k, kmerInterval) {
				if in, _ := bloom.GetBit(uint64(hash)); !in{ 
					srPos += kmerInterval
					continue 
				}
				positions := kmerHash[hash]
				//skip if too many kmers
				if len(positions) < 10 {
					for _, pos := range positions {
						// fmt.Println(srPos, pos)
						// fmt.Println(ref[pos:pos+k])
						aln, offset := getFuzzy(&sr, &ref, srPos, pos, k, window, thres)
						alnLen := len(aln)
						if alnLen > bestAlnLen {
							bestAln, bestAlnLen, bestOffset = aln, alnLen, offset
						}
						if alnLen > 25 {
							prevCandidates = append(prevCandidates, pair{srPos, pos})
						}
					}
				}
				srPos += kmerInterval
			}
		}
		prevSimilar = diff < 0.05
		if bestAlnLen > 15 {
			//do smith-waterman alignment
			printAln(bestAln, bestOffset)
		}
	}
	elapsed := time.Since(start)
    fmt.Println(elapsed)
	
	fmt.Println("done")
}