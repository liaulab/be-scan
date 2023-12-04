package main

import (
  "github.com/TheFutureIsOurs/ahocorasick"
  "bufio"
  "fmt"
  "log"
  "os"
  "math"
  "encoding/json"
  "flag"
  "strconv"
  "strings"
)

func main() {

}

func ref_genome_check(gRNAs_filepath string, genome_filepath string, outputfile string) {

  // build aho corasick automata from gRNAs
  gRNAs, err, gRNA_len := readLines(gRNAs_filepath)
  ac, err := ahocorasick.Build(gRNAs)
  check_err(ac, err)
  // ac, err := ahocorasick.BuildFromFile("./henrykim_gRNAs.txt")

  // open genome file
  file, err := os.Open(genome_filepath)
  if err != nil {
    log.Fatal(err)
    fmt.Printf("Error: %s", err.Error())
  }
  defer file.Close()

  // scan through genome line by line
  scanner := bufio.NewScanner(file)
  const maxCapacity int = 45000000
  buf := make([]byte, maxCapacity)
  scanner.Buffer(buf, maxCapacity)

  // two results are indices of a given gRNA where matches occur, and total count of matches
  var result_indices = make(map[string][]int)
  var result_counts = make(map[string]int)
  var prev_line = strings.Repeat("N", gRNA_len)
  var curr_line = strings.Repeat("N", gRNA_len)
  var ind = 0

  for scanner.Scan() {
    prev_line = curr_line
    curr_line = scanner.Text()
    var slice = math.Min(float64(gRNA_len), float64(len(curr_line)))
    search := ac.MultiPatternSearch([]rune(prev_line+curr_line[:int(slice)]))

    // for every line, scan for gRNAs with aho corasick automata and update maps
    for _, v := range search {
      if _, exists := result_indices[string(v.Value)]; exists {
        result_indices[string(v.Value)] = append(result_indices[string(v.Value)], ind+v.Begin)
        result_counts[string(v.Value)] += 1
      } else {
        result_indices[string(v.Value)] = []int{ind+v.Begin}
        result_counts[string(v.Value)] = 1
      }
    }
    ind += len(curr_line)
  }

  // save maps to json, CHANGE output filenames here
  save_to_json(result_indices, "./" + outputfile + ".json")
  save_to_json(result_counts, "./" + outputfile + "_count.json")

}


func save_to_json(results interface{}, name string) {
  jsonData, err := json.Marshal(results)
  if err != nil {
    fmt.Printf("Error: %s", err.Error())
  }

  jsonFile, err := os.Create(name)
  if err != nil {
    panic(err)
  }
  defer jsonFile.Close()

  jsonFile.Write(jsonData)
  jsonFile.Close()
  fmt.Println("JSON data written to ", jsonFile.Name())
}

func check_err(ac *ahocorasick.Ac, err error) {
  if err != nil {
    log.Fatal(err)
    fmt.Printf("Error: %s", err.Error())
  }
  print(ac)
  print("\n")
}

func readLines(path string) ([]string, error, int) {
    file, err := os.Open(path)
    if err != nil {
        return nil, err, 0
    }
    defer file.Close()

    gRNA_len := 0
    var lines []string
    scanner := bufio.NewScanner(file)
    for scanner.Scan() {
      gRNA := scanner.Text()
      lines = append(lines, gRNA)
      
      gRNA_len = int(math.Max(float64(gRNA_len), float64(len(gRNA))))
    }
    return lines, scanner.Err(), gRNA_len
}
