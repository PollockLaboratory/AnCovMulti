package seqmerAC

import (
	"bufio"
	"fmt"
	"globals"
	"os"
	"strconv"
	"strings"
)

// slices of nucleotide single letter names and complementary names
var nucs = []byte{'G', 'A', 'C', 'T'}
var cnucs = []byte{'C', 'T', 'G', 'A'}

// sequences holds sequence info
type Sequences struct {
	Name      string
	seqfile   string
	seqmap    map[string]string
	seqfilter map[string]*filter
	Filterdir string
	minlength int
	record    bool
	dofilter  bool
	linelimit int
	linemin   int
	count     int
	isalign   bool
	firstlen  int
	totallen  int
	// could add length map
	outfile    string
	filetype   string
	entrystart string
}

// Init creates new parameter structure of hash types
func (seqs *Sequences) Init(seqfile string, minlen int, llimit int, minline int, dorecord bool, filetype string, dofilter bool) {
	seqs.Name = seqfile
	seqs.seqfile = seqfile
	seqs.seqmap = make(map[string]string)
	seqs.seqfilter = make(map[string]*filter)
	seqs.Filterdir = "none"

	seqs.minlength = minlen
	seqs.linelimit = llimit
	seqs.linemin = minline
	seqs.record = dorecord
	seqs.count = 0
	seqs.isalign = false
	seqs.dofilter = dofilter
	fmt.Println("input and recorded filter states", dofilter, seqs.dofilter)
	seqs.firstlen = 0
	seqs.totallen = 0
	seqs.outfile = "sequences.txt"
	seqs.filetype = filetype
	if filetype == "fastq" {
		seqs.entrystart = "@"
	} else {
		seqs.entrystart = ">"
	}
}

// kmers.matches = make([]*oligo, 0)

// filter holds info about read filters
type filter struct {
	name      string
	readlen   int
	pstart    int
	pend      int
	Wuleft    int
	Wuright   int
	fp        int
	lp        int
	direction string

	olist []string // stable ordered list pointing to kmer info
	klist []*oligo // stable ordered list pointing to kmer info
} //

// reads holds info about the reads (but not the sequences)
type reads struct {
	name   string
	length int
	start  int
	end    int
} //

// primer holds start end and direction of primer
type primer struct {
	pstart int
	pend   int
	dir    string
} //

// oligo holds info for individual kmer
type oligo struct {
	name    string
	revcomp string
	poses   []int
	kcount  int
} // will make global seqs

// Init creates new parameter structure of hash types
func (k *oligo) Init(kmer string, kcount int) {
	k.name = kmer
	k.revcomp = rc(kmer)
	k.kcount = kcount
	k.poses = make([]int, 0)
	// fmt.Println("Added k.poses", k.name)
}

// oligos holds kmer map of all kmers, tracks length and total kmers
type Oligos struct {
	name       string
	kmap       map[string]*oligo
	rmap       map[string]*oligo
	kcount     map[string]int
	locprimer  map[int]int
	readcounts map[int]int
	qmatches   map[string]*QSeqMatches
	primers    map[int]*primer
	total      int
	klen       int
	minprint   int
	Outfile    string
	kfile      string
	printNs    bool
	remnant    string
} // will make global kmers

// Init creates new parameter structure of hash types
func (kmers *Oligos) Init(klen int, koutfile string, doNs bool, kminprint int) {
	fmt.Println("Running kinit ", koutfile)
	kmers.name = "generic_seqfile"
	kmers.kmap = make(map[string]*oligo)
	kmers.rmap = make(map[string]*oligo)
	kmers.kcount = make(map[string]int)
	kmers.locprimer = make(map[int]int)
	kmers.readcounts = make(map[int]int)
	kmers.qmatches = make(map[string]*QSeqMatches)
	kmers.primers = make(map[int]*primer)
	kmers.total = 0
	kmers.minprint = kminprint
	kmers.klen = klen
	kmers.Outfile = koutfile
	kmers.printNs = doNs
}

// Clearqmatches turns all qmatches to nil and makes new map
func (kmers *Oligos) Clearqmatches() {
	kmers.qmatches = nil
	kmers.qmatches = make(map[string]*QSeqMatches)
}

// QSeqMatches holds info for matches to reference sequence
type QSeqMatches struct {
	ID       int
	seqname  string
	seqlen   int
	matches  []*oligo
	matchpos []int
} // will be held in sequence info

// # routines

//
// Read and print primers and kmers
//

// Kprint outputs kmer counts
func (kmers *Oligos) Kprint() {
	fmt.Println("Opening Kmer Count Output File", kmers.Outfile)
	fkout, _ := os.Create(kmers.Outfile)
	defer fkout.Close()
	kwriter := bufio.NewWriter(fkout)
	defer kwriter.Flush() // need this to get output

	fmt.Fprintln(kwriter, "kmer\tcount")
	for kmer, kcount := range kmers.kcount { // fix this so works with structure
		if (kcount >= kmers.minprint) && (kmers.printNs || (!strings.Contains(kmer, "N"))) {
			fmt.Fprintf(kwriter, "%s\t%d\n", kmer, kcount)
		}
	}
}

// Kposprint2 outputs kmer positions from reference qmers
//  prints a bit more thean Kposprint()
func (qmers *Oligos) Kposprint2(refmers *Oligos, kmax int) {
	fmt.Println("Opening Kmer Position Output File", qmers.Outfile)
	fkout, _ := os.Create(qmers.Outfile)
	defer fkout.Close()
	kwriter := bufio.NewWriter(fkout)
	defer kwriter.Flush() // need this to get output

	fmt.Fprintln(kwriter, "kmer\tcount\tkmerloc\treadcount\tpositions")
	qmin := qmers.minprint
	for qmer, kcount := range qmers.kcount { // fix this so works with structure
		if (kcount >= qmin) && (qmers.printNs || (!strings.Contains(qmer, "N"))) {
			if refmers.kmap[qmer] != nil {
				poses := refmers.kmap[qmer].poses
				if (poses != nil) && (kcount < kmax) && (kcount >= qmin) {
					kmerloc := poses[0]
					readcount := qmers.readcounts[kmerloc]
					fmt.Fprintf(kwriter, "%s\t%d\t%d\t%d", qmer, kcount, kmerloc, readcount)
					for i := range poses {
						fmt.Fprintf(kwriter, "\t%d", poses[i])
					}
					fmt.Fprintf(kwriter, "\n")
				}
			}
		}
	}
	fmt.Println("\nprinting at the end here jkpq")
}

// Kposprint outputs kmer positions
// used to print the refmers
func (kmers *Oligos) Kposprint(outfile string, kmin int, kmax int) {
	fmt.Println("Opening Kmer Position Output File", outfile)
	fkout, _ := os.Create(outfile)
	defer fkout.Close()
	kwriter := bufio.NewWriter(fkout)
	defer kwriter.Flush() // need this to get output

	fmt.Fprintln(kwriter, "kmer\tcount\tpositions")
	for kmer, kcount := range kmers.kcount { // fix this so works with structure
		if (kcount >= kmin) && (kmers.printNs || (!strings.Contains(kmer, "N"))) {
			if kmers.kmap[kmer] != nil {
				poses := kmers.kmap[kmer].poses
				if (poses != nil) && (kcount < kmax) && (kcount >= kmin) {
					fmt.Fprintf(kwriter, "%s\t%d", kmer, kcount)
					for i := range poses {
						fmt.Fprintf(kwriter, "\t%d", poses[i])
					}
					fmt.Fprintf(kwriter, "\n")
				}
			}
		}
	}
}

// Readk inputs kmer counts
func (kmers *Oligos) Readk(kcountfile string) {
	// reading stuff
	kmers.kfile = kcountfile
	fmt.Println("File to open for Readk() is ", kmers.kfile)
	fpin, err := os.Open(kmers.kfile)
	globals.Check(err)
	defer fpin.Close()
	scanner := bufio.NewScanner(fpin)

	// read, record, count kmers
	const splitter = "\t"
	const bipart = 2
	var linecount int
	for scanner.Scan() {
		if linecount > 0 {
			line := scanner.Text() // should not include eol
			tokens := strings.Split(line, splitter)
			kmer := tokens[0]
			count, _ := strconv.Atoi(tokens[1])
			kmers.kcount[kmer] = count
			kmers.kmap[kmer] = new(oligo)
			kinfo := kmers.kmap[kmer]
			kinfo.name = kmer
			kinfo.kcount = count
			kinfo.revcomp = rc(kmer)
			kmers.rmap[kinfo.revcomp] = kinfo
			kinfo.poses = make([]int, 0) // imagining option to max pos at 10

			//fmt.Println(kmer,"\t",count)
		}
		linecount++
	}
	fmt.Println("line count in Readk is ", linecount)
}

// ReadPrimers gets primer locations
func (kmers *Oligos) ReadPrimers(primerfile string, direction string) {
	fmt.Println("File to open for ReadPrimers() is ", primerfile)
	fpin, err := os.Open(primerfile)
	globals.Check(err)
	defer fpin.Close()
	scanner := bufio.NewScanner(fpin)
	// read, record, primer locations
	const splitter = "\t"
	for scanner.Scan() {
		line := scanner.Text() // should not include eol
		tokens := strings.Split(line, splitter)
		pstart, _ := strconv.Atoi(tokens[0])
		pend, _ := strconv.Atoi(tokens[1])
		// this ought to safety check if already exists
		if direction == "right" {
			temp := pstart
			pstart = pend
			pend = temp
		}
		kmers.primers[pstart] = new(primer)
		primer := kmers.primers[pstart]
		primer.pstart = pstart
		primer.pend = pend
		primer.dir = direction
	}
}

// PrintPrimers prints out list of primers to outfile
func (kmers *Oligos) PrintPrimers(outfile string) {
	fmt.Println("Opening primer Output File", outfile)
	fkout, _ := os.Create(outfile)
	defer fkout.Close()
	kwriter := bufio.NewWriter(fkout)
	defer kwriter.Flush() // need this to get output

	fmt.Fprintln(kwriter, "firstpos\tendpos\tdirection")
	for _, primer := range kmers.primers {
		fmt.Fprintf(kwriter, "%d\t%d\t%s\n", primer.pstart, primer.pend, primer.dir)
	}
}

//
// QueryNotRef functions
//

// QueryNotRef puts qmer counts in qnotkmers if not in kmers
func (qnotkmers *Oligos) QueryNotRef(kmers *Oligos, qmers *Oligos, kqmers *Oligos, direction string) {
	kqcount := 0
	qnkcount := 0
	for qmer, _ := range qmers.kcount {
		origimer := qmer
		if direction == "right" {
			qmer = rc(qmer) // original but rc
			//fmt.Println("flipping directions", origimer, qmer, direction, kmers.kcount[origimer], kmers.kcount[qmer])
		} else {
			//fmt.Println("not flipping", origimer, qmer, direction, kmers.kcount[origimer], kmers.kcount[qmer])
		}
		qmer = strings.ToUpper(qmer)         // counting lower case will only happen when it hits
		if _, ok := kmers.kcount[qmer]; ok { // match to kmers always UC
			kqmers.kcount[origimer] += qmers.kcount[origimer]
			kqcount++
		} else {
			// fmt.Println("putting in qnotk", origimer, qmer, direction, kmers.kcount[origimer], kmers.kcount[qmer])
			qnotkmers.kcount[origimer] += qmers.kcount[origimer]
			qnkcount++
		}
	}
	fmt.Println("total each type", kqcount, qnkcount)
}

//
// Qfilter printing
//

// getpos gets and returns the (first) position for a kmer of interest
func (kmers *Oligos) getpos(kinfo *oligo, kmax int, kmin int, merpos int) int {
	kmer := kinfo.name
	kcount := kmers.kcount[kmer]
	if (kcount >= kmin) && (kmers.printNs || (!strings.Contains(kmer, "N"))) {
		if kmers.kmap[kmer] != nil {
			refposes := kmers.kmap[kmer].poses
			if (refposes != nil) && (kcount < kmax) && (kcount >= kmin) {
				return refposes[0]
			}
		}
	}
	return 0
}

// QFilterPrint outputs info on read matches to reference kmers
//func (kmers *Oligos) QFilterPrint(outfile string, kmin int, kmax int, printheader bool) {
// printheader bool
func (kmers *Oligos) QFilterPrint(globs *globals.Params, outfile string, kmin int, kmax int) {
	fmt.Println("Opening Kmer Position Output File", outfile)
	fkout, _ := os.Create(outfile)
	defer fkout.Close()
	kwriter := bufio.NewWriter(fkout)
	defer kwriter.Flush() // need this to get output

	var firstcounts map[int]int
	firstcounts = make(map[int]int)
	var secondcounts map[int]int
	secondcounts = make(map[int]int)
	var lengthcounts map[int]int
	lengthcounts = make(map[int]int)

	qmseqcount := 0
	goodseqcount := 0
	semigoodseq := 0
	sortagoodseq := 0
	fmt.Fprintln(kwriter, "seqname\tseqlen\tstart\tstop\tfp\tlp\tWuLeft\tWuRight\tdirection\thits\tmisses\tdensity\tWudiff\tposdiff")
	for qname, qmatch := range kmers.qmatches {
		goodseq := true
		qmseqcount++
		klen := kmers.klen
		matchcount := counthits2(qmatch.matchpos)
		if matchcount < globs.Geti("minmatch") {
			goodseq = false
		}
		seqlen := qmatch.seqlen
		if seqlen < globs.Geti("minseq") {
			goodseq = false
		}
		possible := seqlen - klen
		misses := possible - matchcount + 1
		if misses > globs.Geti("maxmisses") {
			goodseq = false
		}
		density := 0.0
		if possible > 0 {
			density = float64(misses) / float64(possible)
		}
		//fmt.Println("Doing qfilter sequence", qname, goodseq, matchcount, misses, seqlen)
		if goodseq {
			firstpos := qmatch.matchpos[0]
			lastpos := qmatch.matchpos[matchcount-1]
			firstmatch := qmatch.matches[0]
			lastmatch := qmatch.matches[matchcount-1]
			orient := "left"
			if firstpos < 0 {
				orient = "right"
			}
			//fmt.Println("Calculated orientation", orient)
			first := kmers.getpos(firstmatch, kmax, kmin, firstpos)
			second := kmers.getpos(lastmatch, kmax, kmin, lastpos)
			if first > second {
				temp := first
				first = second
				second = temp
			}
			firstcounts[first]++
			secondcounts[second]++
			lengthcounts[seqlen]++
			if firstpos < 0 {
				firstpos = -firstpos
				lastpos = -lastpos
			}
			if firstpos > globs.Geti("maxfirsthit") {
				goodseq = false
			}
			if lastpos < globs.Geti("minlasthit") {
				goodseq = false
			}
			posdiff := lastpos - seqlen + kmers.klen - 1
			if posdiff < globs.Geti("minposdiff") {
				goodseq = false
			}
			Wudiff := second - first + klen
			if Wudiff > globs.Geti("Wudiffmax") {
				goodseq = false
			}
			if goodseq {
				sortagoodseq++
			}

			// check on match to primer start
			pstart := 0
			pend := 0
			// works for left
			if orient == "left" {
				primer := kmers.primers[first]
				if primer == nil { // could pull this out of loop
					goodseq = false
				} else {
					pend = primer.pend - first
					pstart = primer.pstart - first + 1
				}
			} else if orient == "right" {
				//fmt.Println("Testing primer", orient)
				primer := kmers.primers[second+klen]
				//fmt.Println("here on primer", orient)
				if primer == nil {
					goodseq = false
				} else {
					pstart = primer.pend - first
					pend = primer.pstart - first
				}
			}
			if goodseq {
				semigoodseq++
			}
			if goodseq {
				goodseqcount++
				fmt.Fprintf(kwriter, "%s\t%d\t%d\t%d\t", qname, seqlen, pstart, pend)
				lastpos = lastpos + kmers.klen - 1
				fmt.Fprintf(kwriter, "%d\t%d\t", firstpos, lastpos)
				fmt.Fprintf(kwriter, "%d\t%d\t", first, second)
				fmt.Fprintf(kwriter, "%s\t%d\t%d\t%f\t", orient, matchcount, misses, density)
				fmt.Fprintf(kwriter, "%d\t%d\t", Wudiff, posdiff)
				fmt.Fprintln(kwriter)
			}
		}
	}
	fmt.Println("First counts\n", "site\tcount")
	fmt.Println("Second counts\n", "site\tcount")
	fmt.Println("\nHow Many?\n", goodseqcount, qmseqcount, sortagoodseq, semigoodseq)
}

//
// Setting up some kmer and seq and filter stuff
//

// GetoutfileMulti creates the right outfile name and adds to kmers
func (kmers *Oligos) GetoutfileMulti(kcountpre string, basename string, directory string) {
	strlen := strconv.Itoa(kmers.klen)
	kmers.Outfile = directory + kcountpre + strlen + "_" + basename + ".xls"
	fmt.Println("AnCov multi outfile name check", kmers.Outfile)
}

// GetSeqnameFilter reads a simple list of seq names on lines for later filtering
func (seqs *Sequences) GetSeqnameFilter(filterfile string) {
	// reading stuff
	fmt.Println("File to open for GetSeqnameFilter() is ", filterfile)
	fpin, err := os.Open(filterfile)
	globals.Check(err)
	defer fpin.Close()
	scanner := bufio.NewScanner(fpin)

	// read sequence names for filter
	const splitter = "\t"
	var linecount int
	direction := seqs.Filterdir
	for scanner.Scan() {
		if linecount >= 0 {
			line := scanner.Text() // should not include eol
			tokens := strings.Split(line, splitter)
			seqname := tokens[0]
			if (len(tokens) > 9) && (tokens[8] == direction) {
				seqs.seqfilter[seqname] = new(filter)
				filter := seqs.seqfilter[seqname]
				filter.readlen, _ = strconv.Atoi(tokens[1])
				filter.pstart, _ = strconv.Atoi(tokens[2])
				filter.pend, _ = strconv.Atoi(tokens[3])
				filter.fp, _ = strconv.Atoi(tokens[4])
				filter.lp, _ = strconv.Atoi(tokens[5])
				filter.Wuleft, _ = strconv.Atoi(tokens[6])
				filter.Wuright, _ = strconv.Atoi(tokens[7])
				filter.direction = tokens[8]
			}
		}
		linecount++
	}
}

// Setfilterdir changes the filter direction for seqs
func (seqs *Sequences) Setfilterdir(direction string) {
	seqs.Filterdir = direction
}

// SetReadCounts adds read counts for every position
func (kmers *Oligos) SetReadCounts(outfile string, direction string) {
	fmt.Println("Opening read counts Output File", outfile)
	fkout, _ := os.Create(outfile)
	defer fkout.Close()
	kwriter := bufio.NewWriter(fkout)
	defer kwriter.Flush() // need this to get output

	movement := 1
	seqlenmax := 255
	if direction == "right" {
		movement = -1
	}
	fmt.Println("Printing kmer location counts to ", outfile)
	fmt.Fprintf(kwriter, "ploc\tpcount\n")
	for ploc, pcount := range kmers.locprimer {
		for i := 0; i < seqlenmax; i++ {
			kmerloc := ploc + (i * movement)
			kmers.readcounts[kmerloc] += pcount
		}
		//fmt.Fprintf(kwriter, "\n%d\t%d\n", ploc, pcount)
		for i := 0; i < seqlenmax; i++ {
			kmerloc := ploc + (i * movement)
			fmt.Fprintf(kwriter, "%d\t%d\n", kmerloc, kmers.readcounts[kmerloc])
		}
	}
}

// MoveReadCounts moves read counts from frommers to kmers
func (kmers *Oligos) MoveReadCounts(frommers *Oligos, outfile string) {
	fmt.Println("Opening read counts Output File", outfile)
	fkout, _ := os.Create(outfile)
	defer fkout.Close()
	kwriter := bufio.NewWriter(fkout)
	defer kwriter.Flush() // need this to get output

	fmt.Println("Printing kmer location counts to ", outfile)
	fmt.Fprintf(kwriter, "ploc\tpcount\n")
	for kmerloc, readcount := range frommers.readcounts {
		kmers.readcounts[kmerloc] = readcount
		fmt.Fprintf(kwriter, "%d\t%d\n", kmerloc, kmers.readcounts[kmerloc])
	}
}

//
// Three sequence readers that do slightly different things
// Tardigrader, Tactopoda, Onychophora
//

//
// Tardigrader functions
//

// Countmers counts kmers in string, adds to stored counts in kmers
// called by Tardigrader()
func (kmers *Oligos) Countmers(seqs *Sequences, seq string, recording bool, name string) {
	if !recording {
		var kmer string
		var filtering bool
		if seqs.seqfilter[name] == nil {
			filtering = false
		} else {
			filtering = seqs.dofilter
		}
		startread := 0
		endread := 10000 // I suppose this could be a problem for v long reads
		Wuleft := 0
		Wuright := 0
		//lp := 0
		fp := 1
		filterdirection := seqs.Filterdir
		direction := seqs.Filterdir
		if filtering {
			direction = seqs.seqfilter[name].direction
			Wuleft = seqs.seqfilter[name].Wuleft
			Wuright = seqs.seqfilter[name].Wuright
			startread = seqs.seqfilter[name].pend // start after end of primer
			if direction == "right" {
				startread = seqs.seqfilter[name].pend - seqs.seqfilter[name].pstart + 1
				//startread = 0
				//endread = seqs.seqfilter[name].pstart
			}
			fp = seqs.seqfilter[name].fp
			if filterdirection == "left" {
				kmers.locprimer[Wuleft]++
			}
			if filterdirection == "right" {
				kmers.locprimer[Wuright]++
			}
		}
		// fmt.Println("in countmers ", seqs.Name, filtering, direction)
		for i := 0; i < (len(seq) - kmers.klen + 1); i++ {
			if (i >= startread) && (i < endread) && (direction == filterdirection) {
				kmer = seq[i : i+kmers.klen]
				if kmers.kmap[kmer] == nil {
					kmers.kmap[kmer] = new(oligo)
					kmers.kmap[kmer].Init(kmer, kmers.kcount[kmer])
				}
				kinfo := kmers.kmap[kmer]
				if filterdirection == "left" {
					kinfo.poses = append(kinfo.poses, (Wuleft + i + fp - 1))
				}
				if filterdirection == "right" {
					kinfo.poses = append(kinfo.poses, (Wuright - i + fp - 1))
				}
				kmers.kcount[kmer]++
				kmers.total++
			} else {
				//println("failed read", i, startread, endread, direction, filterdirection)
			}
		}
		nextpos := len(seq) - kmers.klen + 1
		kmers.remnant = seq[nextpos : nextpos+kmers.klen-1]
	}
}

// Tardigrader reads fasta or fastq file and records into seqs structure or counts kmers
// this is a trimmed down version to do as little as possible
// additionally, this makes use of the filters
func (seqs *Sequences) Tardigrader(kmers *Oligos) {
	var name, seq string
	var count, lcount int

	// reading stuff
	fmt.Println("File to open is thyq2", seqs.seqfile)
	fpin, err := os.Open(seqs.seqfile)
	globals.Check(err)
	defer fpin.Close()
	scanner := bufio.NewScanner(fpin)
	entrycount := 0 // track line in each entry
	entrylimit := 0 // default every line in entry is counted
	// hackish fastq, just reading the first line after name line
	if seqs.filetype == "fastq" {
		entrylimit = 1
	}
	fmt.Println("File type is ", seqs.filetype, "and entry limit is", entrylimit)
	//	fmt.Println("dofilter ", seqs.dofilter, "and filterfile", seqs.seqfilter)
	fmt.Println("In Targdigrader, dofilter ", seqs.dofilter)
	passfilter := false // flag to see if name is in filter list

	// read, record, count kmers
	for scanner.Scan() {
		lcount += 1
		line := scanner.Text()              // should not include eol
		trimline := strings.TrimSpace(line) // trim off leading and lagging whitespace
		nofilter := !seqs.dofilter
		if strings.HasPrefix(line, seqs.entrystart) {
			// fmt.Println("header line", line)
			name = strings.TrimPrefix(trimline, seqs.entrystart)
			count += 1
			entrycount = 1
			kmers.remnant = ""
			if (seqs.seqfilter[name] != nil) && seqs.dofilter {
				//fmt.Println("New seq", name, "number", count)
				passfilter = true
				//fmt.Println("criteria in names", nofilter, passfilter, lcount, entrylimit, entrycount)
			} else {
				passfilter = false
				if (count % 10000) == 0 {
					fmt.Println("Skip seq", name, "number", count)
				}
			}
		} else {
			//fmt.Println("criteria out of names", nofilter, passfilter, lcount, seqs.linelimit, entrylimit, entrycount)
			if passfilter {
				//fmt.Println("Testing seq ", name, "number", count)
				//fmt.Println("criteria ", nofilter, passfilter, lcount, entrylimit, entrycount)
			}
			if (nofilter || passfilter) && (lcount < seqs.linelimit) && (entrylimit < 1 || entrycount <= entrylimit) {
				//fmt.Println("Doing seq ", name, "number", count)
				seq = kmers.remnant + trimline
				entrycount++
				if (lcount <= seqs.linelimit) && (lcount > seqs.linemin) {
					kmers.Countmers(seqs, seq, seqs.record, name)
				}
				if lcount < 50000 {
					if (lcount % 10000) == 0 {
						fmt.Println(lcount, len(seq))
					}
				}
				if (lcount % 50000) == 0 {
					fmt.Println(lcount, len(seq))
				}
			}
		}
	}
	fmt.Println("Lines counted\n", count, lcount)
}

//
// Tactopoda functions
//

// addtoklist adds pointer to existing kinfo
// called by locater()
func (seqs *Sequences) addtoklist(seq string, pos int, kinfo *oligo, orient string) {
	if seqs.seqfilter[seq] == nil {
		seqs.seqfilter[seq] = new(filter)
	}
	if seqs.seqfilter[seq].klist == nil {
		seqs.seqfilter[seq].klist = make([]*oligo, 0)
	}
	if seqs.seqfilter[seq].olist == nil {
		seqs.seqfilter[seq].olist = make([]string, 0)
	}
	filter := seqs.seqfilter[seq]
	filter.olist = append(filter.olist, orient)
	filter.klist = append(filter.klist, kinfo)
}

// Locater locates kmers in string, adds to stored counts in kmers
// called by Tactopoda()
func (kmers *Oligos) Locater(seqs *Sequences, seq string, seqname string, recording bool, linepos int, kmax int, kmin int) {
	if !recording {
		var token string
		for i := 0; i < (len(seq) - kmers.klen + 1); i++ {
			token = seq[i : i+kmers.klen]
			kinfo := kmers.kmap[token] // have to reassign because out of loop
			if kinfo != nil {
				if (kinfo.kcount < kmax) && (kinfo.kcount >= kmin) {
					kinfo.poses = append(kinfo.poses, (linepos + i))
				}
				seqs.addtoklist(seqname, i, kinfo, "left")
			} // dies when there is an N, so make sure not nil pointer
		}
		nextpos := len(seq) - kmers.klen + 1
		kmers.remnant = seq[nextpos : nextpos+kmers.klen-1]
	}
}

// Tactopoda reads fasta or fastq file and records positions of kmers
// note: this was previously passing the neighbors structure but only using the max
func (seqs *Sequences) Tactopoda(kmers *Oligos, kmax int, kmin int, outfile string) {
	var name string
	var seq string
	var count, lcount int

	// reading stuff
	fmt.Println("File to open is ", seqs.seqfile)
	fpin, err := os.Open(seqs.seqfile)
	globals.Check(err)
	defer fpin.Close()
	scanner := bufio.NewScanner(fpin)

	// writing positions
	fkout, _ := os.Create(outfile)
	defer fkout.Close()
	kwriter := bufio.NewWriter(fkout)
	defer kwriter.Flush() // need this to get output
	fmt.Fprintf(kwriter, "%s\t%s\t%s\n", "New seq", "number", "position")

	entrycount := 0 // track line in each entry
	entrylimit := 0 // default every line in entry is counted
	position := 1   // starting position is 1
	// hackish fastq, just reading the first line after name line
	if seqs.filetype == "fastq" {
		entrylimit = 1
	}
	fmt.Println("File type is ", seqs.filetype, "and entry limit is", entrylimit)
	// read, record, count kmers
	for scanner.Scan() {
		lcount += 1
		line := scanner.Text()              // should not include eol
		trimline := strings.TrimSpace(line) // trim off leading and lagging whitespace
		if strings.HasPrefix(line, seqs.entrystart) {
			fmt.Println("header line", line)
			name = strings.TrimPrefix(trimline, seqs.entrystart)
			count += 1
			entrycount = 1
			kmers.remnant = ""
			if (lcount <= seqs.linelimit) && (lcount > seqs.linemin) {
				//fmt.Fprintf(kwriter, "%s\t%d\t%d\n", name, count, (position - kmers.klen))
			}
		} else if entrylimit < 1 || entrycount <= entrylimit {
			seq = kmers.remnant + trimline // this is why correct seq len in pos
			entrycount++
			if (lcount <= seqs.linelimit) && (lcount > seqs.linemin) {
				kmers.Locater(seqs, seq, name, seqs.record, position, kmax, kmin)
				position += (len(seq) - kmers.klen + 1)
				//fmt.Println("updating position",seq, len(seq),name,"and running locater at",position)
			}
			if lcount < 100000 {
				if (lcount % 10000) == 0 {
					fmt.Println(lcount, len(seq))
				}
			}
			if (lcount % 100000) == 0 {
				fmt.Println(lcount, len(seq))
			}
		}
	}
	fmt.Println("Lines counted\n", count, lcount)
}

//
// Onychophora functions
//

// addmatch adds the match information to kmers qmatches struct
// if no hit in kmap, looks for match in rmap
// returns integer number of hits found
// called by matchmer
func (kmers *Oligos) addmatch(kmer string, qseqname string, readpos int) int {
	qmatch := kmers.qmatches[qseqname]
	kinfo := kmers.kmap[kmer]
	var nohits int
	if kinfo != nil {
		qmatch.matches = append(qmatch.matches, kinfo)
		qmatch.matchpos = append(qmatch.matchpos, readpos)
	} else {
		kinfo = kmers.rmap[kmer]
		if kinfo != nil {
			qmatch.matches = append(qmatch.matches, kinfo)
			qmatch.matchpos = append(qmatch.matchpos, -readpos)
		}
	}
	if kinfo != nil {
		nohits = 0
	} else {
		nohits = 1
	}
	return nohits
}

// addqmatch adds the qmatch information for a sequence to kmers qmatches struct
// called by matchmer
func (kmers *Oligos) addqmatch(name string, seqnum int, seqlen int) {
	if _, ok := kmers.qmatches[name]; !ok {
		kmers.qmatches[name] = new(QSeqMatches)
		kmers.qmatches[name].matches = make([]*oligo, 0)
		kmers.qmatches[name].matchpos = make([]int, 0)
	}
	qmatch := kmers.qmatches[name]
	qmatch.ID = seqnum
	qmatch.seqname = name
	qmatch.seqlen = seqlen
} // this is still lazily redoing the ID and seqname

// matchmer find kmers in string that match stored kmer file
// called by Onychophora
func (kmers *Oligos) matchmer(seq string, name string, seqnum int) {
	var token string
	seqlen := len(seq)
	kmers.addqmatch(name, seqnum, seqlen)
	nohits := 0
	for i := 0; i < (seqlen - kmers.klen + 1); i++ {
		token = seq[i : i+kmers.klen]
		nohits += kmers.addmatch(token, name, i+1)
	}
	if nohits > 0 {
		for i := 0; i < (seqlen - kmers.klen + 1); i++ {
			token = seq[i : i+kmers.klen]
		}
	}
	nextpos := len(seq) - kmers.klen + 1
	kmers.remnant = seq[nextpos : nextpos+kmers.klen-1]
}

// Onychophora (velvet worms) reads fasta or fastq file and records positions of kmers
// note: this was previously passing the neighbors structure but only using the max
// calls matchmer()
func (seqs *Sequences) Onychophora(kmers *Oligos, kmax int, kmin int, outfile string) {
	var name, seq string
	var count, lcount int

	// reading stuff
	fmt.Println("File to open is ", seqs.seqfile, seqs.filetype)
	fpin, err := os.Open(seqs.seqfile)
	globals.Check(err)
	defer fpin.Close()
	scanner := bufio.NewScanner(fpin)

	// writing positions
	fkout, _ := os.Create(outfile)
	defer fkout.Close()
	kwriter := bufio.NewWriter(fkout)
	defer kwriter.Flush() // need this to get output
	fmt.Fprintf(kwriter, "%s\t%s\t%s\n", "New seq", "number", "postion")

	entrycount := 0 // track line in each entry
	entrylimit := 0 // default every line in entry is counted
	position := 1   // starting position is 1
	// hackish fastq, just reading the first line after name line
	if seqs.filetype == "fastq" {
		entrylimit = 1
	}
	fmt.Println("File type is ", seqs.filetype, "and entry limit is", entrylimit)
	// read, record, count kmers
	for scanner.Scan() {
		lcount += 1
		line := scanner.Text()              // should not include eol
		trimline := strings.TrimSpace(line) // trim off leading and lagging whitespace
		if strings.HasPrefix(line, seqs.entrystart) {
			// fmt.Println("header line", line)
			name = strings.TrimPrefix(trimline, seqs.entrystart)
			count += 1
			entrycount = 1
			kmers.remnant = ""
			// fmt.Println("more info", name, count)
			if (lcount <= seqs.linelimit) && (lcount > seqs.linemin) {
				// fmt.Fprintf(kwriter, "%s\t%d\t%d\n", name, count, (position - kmers.klen))
			}
		} else if entrylimit < 1 || entrycount <= entrylimit {
			seq = kmers.remnant + trimline // this is why correct seq len in pos
			entrycount++
			if (lcount <= seqs.linelimit) && (lcount > seqs.linemin) {
				kmers.matchmer(seq, name, count)
				position += (len(seq) - kmers.klen + 1)
				//fmt.Println("updating position",seq, len(seq),name,"and running locater at",position)
			}
			keepcompany(lcount, len(seq))
		}
	}
	fmt.Println("Lines counted\n", count, lcount)
}

//
// sequence manipulation basics
//

// rc returns reverse complement
func rc(kmer string) string {
	kbits := []byte(kmer)
	rcbits := []byte(kmer)
	klen := len(kmer)
	for i := range kbits {
		for n := range nucs {
			if nucs[n] == kbits[i] {
				rcbits[klen-i-1] = cnucs[n]
			}
		}
	}
	return string(rcbits)
}

// recordseq adds seq, name, and count total information to seqs
func (seqs *Sequences) recordseq(seq string, name string) {
	seqlen := len(seq)
	if (seqlen > seqs.minlength) && seqs.record {
		blurb := "Error 73, sequence name already exists => "
		if _, exists := seqs.seqmap[name]; exists {
			panic(blurb + name)
		}
		// a fancier version would add a number to the name or something instead of panicking
		seqs.seqmap[name] = seq
		seqs.count += 1
		seqs.totallen += len(seq)
		fmt.Println("Recorded sequence ", name, "length", seqlen)
	}
}

// stats basics

// counthits2 returns the total non-zero numbers in an integer slice
func counthits2(intslice []int) int {
	var hits int
	if intslice != nil {
		for i := range intslice {
			if intslice[i] != 0 {
				hits += 1
			}
		}
	}
	return hits
}

// counthits returns the total positive numbers in an integer slice
func counthits(intslice []int) int {
	var hits int
	if intslice != nil {
		for i := range intslice {
			if intslice[i] > 0 {
				hits += 1
			}
		}
	}
	return hits
}

// getsum returns the sum of elements in an integer slice
func getsum(intslice []int) int {
	var sum int
	if intslice != nil {
		for i := range intslice {
			if intslice[i] > 0 {
				sum += intslice[i]
			}
		}
	}
	return (sum)
}

// sumdiffsqr returns the sum of squared differences of the mean of integer slice
func sumdiffsqr(intslice []int, mean float64) float64 {
	var sumsqr, realval, diff float64
	if intslice != nil {
		for i := range intslice {
			if intslice[i] > 0 {
				realval = float64(intslice[i])
				diff = realval - mean
				sumsqr += diff * diff
			}
		}
	}
	return (sumsqr)
}

// keepcompany outputs count and seqlen at intervals
func keepcompany(count int, seqlen int) {
	if count < 100000 {
		if (count % 10000) == 0 {
			fmt.Println(count, seqlen)
		}
	}
	if (count % 100000) == 0 {
		fmt.Println(count, seqlen)
	}
}
