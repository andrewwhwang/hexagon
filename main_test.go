package main

import (
	"fmt"
	"sync"
	"testing"
)

var ref = "CACATGGAATTAGGCCAGTAGTATCAACTCAACTGCTGTTAAATGGCAGTCTAGCAGAAGAAGAGGTAGTAATTAGATCTGCCAATTTCACAGACAATGCTAAAACCATAATAGTACAGCTGAACCAATCTGTAGAAATTAATTGTACAAGACCCAACAACAATACAAGAAAAAGAATCCGTATCCAGAGGGGACCAGGGAGAGCATTTGTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGTAACATTAGTAGAGCAAAATGGAATGCCACTTTAAAACAGATAGCTAGCAAATTAAGAGAACAATTTGGAAATAATAAAACAATAATCTTTAAGCAATCCTCAGGAGGGGACCCAGAAATTGTAACGCACAGTTTTAATTGTGGAGGGGAATTTTTCTACTGTAATTCAACACAACTGTTTAATAGTACTTGGTTTAATAGTACTTGGAGTACTGAAGGGTCAAATAACACTGAAGGAAGTGACACAATCACACTCCCATGCAGAATAAAACAATTTATAAACATGTGGCAGGAAGTAGGAAAAGCAATGTATGCCCCTCCCATCAGCGGACAAATTAGATGTTCATCAAATATTACAGGGCTGCTATTAACAAGAGATGGTGGTAATAACAACAATGAGTCCGAGATCTTCAGACCTGGAGGAGGAGATATGAGGGACAATTGG"
var sr = "AAAAGAACAACTTGGAAATAAAACAATAGTCTTTAATCAATCCTCAGGAGGGGACCCAGAAATTGTACTGCACAATTTTAATTGTAGAGGGGAATTTTTCTACTGTAATCTATCGCAACTGTTTAATTGGACAGATAATGAGGCTTGGAGTAATGAAACCTATACTACTGACACTGTGATAGTGTCAGTAGTATTGGTTTCATTACTCCAAGCCTCATTATCTGTCCAATTAAACAGTTGCGATAGATTACAGTAGAAAAATTCCCCTCTACAATTAAAATTGTGCAGTACAATTTCTGGGTCCCCTCCTGAGGATTGATTAAAGACTATTGTTTTATTTCCAAGTTGTT"

func TestGetFuzzy(t *testing.T) {
	// fmt.Println("asdfasdf")
}

func TestExtendHead(t *testing.T) {
	var wg sync.WaitGroup
	var start int
	wg.Add(1)
	srPos := 184
	refPos := 14
	k := 8
	window := 5
	thres := 3

	go func() {
		start = extendHead(sr, &ref, srPos, refPos, k, thres, window)
		wg.Done()
	}()
	wg.Wait()
	fmt.Println(start)
}

// func TestExtendTail(t *testing.T) {
// 	var wg sync.WaitGroup
// 	var end int
// 	wg.Add(1)
// 	srPos := 28
// 	refPos := 215
// 	k := 8
// 	window := 8
// 	thres := 5

// 	go func() {
// 		end = extendTail(sr, &ref, srPos, refPos, k, thres, window)
// 	}()
// 	wg.Wait()
// 	fmt.Println(end)
// }

func TestReverse(t *testing.T) {
	reverse := reverse("abcdefg")
	// t.Log(reverse)
	fmt.Println(reverse)
}
