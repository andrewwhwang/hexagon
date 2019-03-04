package main

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"github.com/ljfuyuan/suffixtree"
)

func getReads(filepath string) chan string{

	file, err := os.Open(filepath)
    if err != nil {
        log.Fatal(err)
    }
	
	ch := make(chan string)
	scanner := bufio.NewScanner(file)

	go func() {
        defer close(ch)
		defer file.Close()
		for i := 0; scanner.Scan(); i++ {
			if i % 4 == 1 {
				// fmt.Println(scanner.Text())
				ch <- scanner.Text()
			}
		}
	}()

	return ch
}

func main() {
	tree := suffixtree.NewGeneralizedSuffixTree()
	index := 0
	for word := range getReads("banana.txt") {
		tree.Put(word, index)
		index++
	}
	ch := make(chan string)

	go func() {
        defer close(ch)
		tree.GetReadsCh(tree.GetRoot(), "", &ch)
	}()
	
	for word := range ch {
		fmt.Println(word)
	}
	fmt.Println("done")
}


// func main() {
// 	words := []string{"banana", "panama"}
// 	tree := suffixtree.NewGeneralizedSuffixTree()
// 	for k, word := range words {
// 		tree.Put(word, k)
// 	}
// 	// indexs := tree.Search("a", -1)

// 	// fmt.Println(indexs)
// 	// //[0 2 1]
// 	// for _, index := range indexs {
// 	// 	fmt.Println(words[index])
// 	// }
// 	suffixtree.PrintNode("", tree.GetRoot())
// 	fmt.Println("done")
// }