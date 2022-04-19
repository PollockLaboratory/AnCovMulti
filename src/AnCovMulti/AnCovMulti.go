package main

import (
	"bufio"
	"fmt"
	"globals"
	"io"
	"os"
	"path"
	"path/filepath"
	seqmer "seqmerAC"
	"strings"
	"time"
)

// program holds information about this program
type program struct {
	id       int
	name     string
	nickname string
	authors  string
	began    string
	modified string
	uses     string
	runrec   string
	computer string
}

var prog program // need a variable to hold the program structure

// these are the files that hold parameters and control the program settings
const factfile = "factory"
const controlfile = "control"
const modefile = "mode"

// Set adds and prints program information
func (p *program) Set(writer io.Writer) {
	p.name = "AnCovMulti" // read in kmer file, check sequences for kmers in seqs
	// optional only consider regional file; report in or not in kmer file; counts from seqs file
	fmt.Fprintln(writer, "\n\tRunning Program", p.name)
	p.authors = "David Pollock"
	p.began = "November 11, 2020"
	p.modified = "December 28, 2020"
	fmt.Fprintln(writer, "\tAuthors:", p.authors, "Last Modified:", p.modified, "\n")
}

func main() {
	fmt.Println("Setting Append File")
	fappend, _ := os.OpenFile("access.log", os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0644)
	defer fappend.Close()
	writer := bufio.NewWriter(fappend)
	defer writer.Flush() // need this to get output

	// set up the program and globals by reading controls
	prog.Set(writer)
	var globs = globals.New()
	globs.ProgSetUp(controlfile, modefile) // should change through factory and command line
	fmt.Fprintln(writer, "The program was run on", globs.Runstart)
	globs.Print(os.Stdout, "\nStatus after Setup")

	// cycle through input files

	// get set of query sequence file names
	qseqroute := getallext(globs.Gets("seqext"), globs.Gets("seqdir"))
	qseqfiles, err := filepath.Glob(qseqroute)
	printerr(err)
	fmt.Println("files that match route", qseqroute, qseqfiles)
	printmatches(qseqfiles)

	// get reference file information
	refdir := globs.Gets("refdir")
	refseqfile := refdir + globs.Gets("seqfile")
	kinfile := refdir + globs.Gets("kinfile")
	left := globs.Gets("left")
	right := globs.Gets("right")
	primeleftfile := refdir + left + globs.Gets("primeroutend")
	primerightfile := refdir + right + globs.Gets("primeroutend")
	qfilterdir := globs.Gets("filterdir")
	outdir := globs.Gets("outputdir")

	fmt.Println("ref seq:\t", refseqfile)
	fmt.Println("kcounts infile:\t", kinfile)
	fmt.Println("left primer bounds:\t", primeleftfile)
	fmt.Println("right primer bounds:\t", primerightfile)
	fmt.Println("q filter directory:\t", qfilterdir)

	// get some kmer variables, mostly about printing
	klen := globs.Geti("klen")
	kcountprefix := globs.Getf("kcountfile")
	kcountinfile := refdir + kcountprefix
	printNs := globs.Getb("printNs")
	kiminprint := globs.Geti("kminprint")
	qminprint := globs.Geti("qminprint")

	// seq file reading variables
	minseqlen := globs.Geti("minseqlen")
	linelimit := globs.Geti("linelimit")
	minline := globs.Geti("minline")
	recseq := globs.Getb("recordseq")
	reftype := globs.Gets("filetype")
	dofilter := globs.Getb("dofilter")

	// these are a bit not necessary but don't want to stop and change
	kmaxprint := globs.Geti("kmaxprint")
	kminprint := globs.Geti("kminprint")
	chromfile := globs.Getf("chrommapfile")
	kposfile := globs.Getf("kpositionsfile")

	// get query location and output information
	qtype := globs.Gets("qfiletype")
	seqnamefilter := refdir + globs.Getf("seqnamefilter") // really need to delete this here, not using it
	kcountdir := globs.Gets("kcountdir")
	pcountdir := globs.Gets("pcountdir")
	prefiltername := globs.Gets("prefilter")
	filterdir := globs.Gets("filterdir")

	// get qnotk output location and file information
	qnotkdir := globs.Getf("qnotkdir")
	kandqdir := globs.Getf("kandqdir")
	qnotkpre := globs.Getf("qnotk")
	kandqpre := globs.Getf("kandq")
	runorient := globs.Gets("runorient")
	bigoutput := globs.Getb("bigoutput")

	//
	// main program here
	//

	fmt.Println("\nStarting main program. qqz\n")

	// read reference kmers and primers, print to validate
	refmers := new(seqmer.Oligos) // new kmers for reference
	refmers.Init(klen, kcountinfile, printNs, kiminprint)
	refmers.Readk(kinfile) // read kmer and counts
	refmers.ReadPrimers(primeleftfile, left)
	refmers.ReadPrimers(primerightfile, right)
	refmers.PrintPrimers(outdir + "primers_out.xls")

	// create seqs and read in reference seqfile
	seqs := new(seqmer.Sequences) // create global; should be a pointer
	seqs.Init(refseqfile, minseqlen, linelimit, minline, recseq, reftype, dofilter)
	seqs.Tactopoda(refmers, kmaxprint, kminprint, chromfile)
	refmers.Kposprint(outdir+kposfile, kminprint, kmaxprint)

	// cycle through query sequence files
	fmt.Println("\nStarting cycling query sequences. qqw\n")
	if globs.Getb("qnotkmode") != true {
		for count, seqfilename := range qseqfiles {
			_, file := path.Split(seqfilename)
			firstsplit := strings.SplitN(file, ".", 2)
			fileID := firstsplit[0]
			if count < globs.Geti("limitseqs") {
				// this first part is just setting up names and outputs
				fmt.Println("\nmatch ", fileID, file)
				query := new(seqmer.Sequences) // create global; should be a pointer
				query.Init(seqfilename, minseqlen, linelimit, minline, recseq, qtype, dofilter)
				query.GetSeqnameFilter(seqnamefilter) // read file with subset of names to analyze

				qmers := new(seqmer.Oligos) // create global; should be a pointer
				qmers.Init(klen, kcountprefix, printNs, kiminprint)
				qmers.GetoutfileMulti(kcountprefix, fileID, outdir+kcountdir) // set up outfile name
				qfiltername := getQFname(filterdir, prefiltername, fileID, ".xls")
				// end setup

				// add query matches to refmers, print filer, then clear matches
				query.Onychophora(refmers, kmaxprint, kminprint, chromfile)
				refmers.QFilterPrint(globs, qfiltername, kminprint, kmaxprint)
				refmers.Clearqmatches()
			}
			fmt.Println("\nRan x out of y query sequences\n", count, globs.Geti("limitseqs"))
		}
	}

	// get set of query sequence file names
	filterroute := getallext("xls", filterdir)
	filterfiles, err := filepath.Glob(filterroute)
	printerr(err)
	fmt.Println("files that match filter route", filterroute)
	// clear and read in refmers again to be safe (used by QueryNotRef)
	refmers = nil
	refmers = new(seqmer.Oligos) // new kmers for reference
	refmers.Init(klen, kcountinfile, printNs, kiminprint)
	refmers.Readk(kinfile) // read kmer and counts

	// cycle through query sequence files along with filterfiles
	fmt.Println("\n Cycle through query sequence files and filterfiles\n")
	for count, seqfilename := range qseqfiles {
		_, file := path.Split(seqfilename)
		firstsplit := strings.SplitN(file, ".", 2)
		fileID := firstsplit[0]
		if count < globs.Geti("limitseqs") {
			fmt.Println("\nmatch ", fileID, file)
			filterfilename := getmatch(fileID, filterfiles, prefiltername)
			if checkFileExists(filterfilename) {
				// set up query sequence and correct filter, create qmers
				query := new(seqmer.Sequences) // create global; should be a pointer
				doqueryfilter := true
				query.Init(seqfilename, minseqlen, linelimit, minline, recseq, qtype, doqueryfilter)
				query.Setfilterdir(runorient)
				query.GetSeqnameFilter(filterfilename) // read file with filtering information

				qmers := new(seqmer.Oligos) // create global; should be a pointer
				qmers.Init(klen, kcountprefix, printNs, kiminprint)
				fmt.Println("\nFilter direction", query.Filterdir)
				qmers.GetoutfileMulti(query.Filterdir+"_"+kcountprefix, fileID, outdir+kcountdir) // set up outfile name
				fmt.Println("qmers outfile name check", qmers.Outfile, query.Name)
				qnkmers := new(seqmer.Oligos) // create global; should be a pointer
				qnkmers.Init(klen, kcountprefix, printNs, qminprint)
				kqmers := new(seqmer.Oligos) // create global; should be a pointer
				kqmers.Init(klen, kcountprefix, printNs, qminprint)

				query.Tardigrader(qmers) // this is the trimmed down count reader
				qnkmers.QueryNotRef(refmers, qmers, kqmers, query.Filterdir)
				qmers.SetReadCounts(outdir+pcountdir+query.Filterdir+"_"+fileID+"pcount.xls", query.Filterdir)
				qnkmers.GetoutfileMulti(query.Filterdir+"_"+qnotkpre, fileID, outdir+qnotkdir) // set up outfile name
				qnkmers.MoveReadCounts(qmers, outdir+pcountdir+query.Filterdir+"_"+fileID+qnotkpre+"pcount.xls")
				qnkmers.Kposprint2(qmers, kmaxprint)

				if bigoutput {
					fmt.Println("Printing bit output files.", bigoutput)
					qmers.Kposprint2(qmers, kmaxprint)

					kqmers.MoveReadCounts(qmers, outdir+pcountdir+query.Filterdir+"_"+fileID+kandqpre+"pcount.xls")
					kqmers.GetoutfileMulti(query.Filterdir+"_"+kandqpre, fileID, outdir+kandqdir)
					kqmers.Kposprint2(qmers, kmaxprint)
				}

			} else {
				fmt.Println("\nNo filter file available for ", fileID)
			}
		}
	}

	// end main code

	fmt.Println("\n\nMain program completed. qqx\n")
	globs.Print(writer, "\nStatus after Program Completion")
	now := time.Now()
	globs.Delta()
	fmt.Fprintln(writer, "\nThe program finished", now)
}

// checkFileExists checks if file exists
func checkFileExists(filename string) bool {
	if filename != "" {
		fmt.Println("We have a file", filename)
		return true
	} else {
		fmt.Println("File is missing", filename)
		return false
	}
}

// getQFname gets name and prints test
func getQFname(dir string, prename string, ID string, append string) string {
	name := dir + prename + ID + append
	fmt.Println("qfilter output name check", name)
	return name
}

func getmatch(infileID string, filenames []string, qfilterpre string) string {
	returnname := "nope"
	for _, filename := range filenames {
		//fmt.Println("testing", infileID, filename)
		_, file := path.Split(filename)
		firstsplit := strings.SplitN(file, ".", 2)
		prelimID := firstsplit[0]
		secondsplit := strings.SplitN(prelimID, qfilterpre, 2)
		outfileID := secondsplit[1]
		//fmt.Println("get match", infileID, outfileID, prelimID)
		if infileID == outfileID {
			//fmt.Println("Match Identified", filename)
			return filename
		}
	}
	return returnname
}

func printmatches(matches []string) {
	for count, s := range matches {
		dir, file := path.Split(s)
		//fmt.Println("match ", count+1, dir, file)
		firstsplit := strings.SplitN(file, ".", 2)
		baserun := firstsplit[0]
		extend := firstsplit[1]
		// fmt.Println("\tfirst split", baserun, extend)
		secondsplit := strings.SplitN(baserun, "_", 3)
		runID := secondsplit[0]
		ID := secondsplit[1]
		direction := secondsplit[2]
		base := runID + "_" + ID
		//fmt.Println("\tsecond split", runID, ID, direction)
		fmt.Println(count+1, base, direction, extend, dir, file, s)
		//fmt.Println("\tfile details", dir, file, s)
	}
}

func getallext(extension string, directory string) string {
	//extension := "fq"
	localdir := "./"
	//directory := "seqfiles/"
	dirpath := localdir + directory
	allfiles := "*"
	pattern := allfiles + "." + extension
	route := dirpath + pattern
	return route
}

func printerr(err error) {
	if err != nil {
		fmt.Println("File reading error", err)
	}
}
