This tutorial describes the functionality of [AnCovMulti v0.1](https://github.com/PollockLaboratory/AnCovMulti)
# Summary
[AnCovMulti](https://github.com/PollockLaboratory/AnCovMulti) reads in multiple sequence read files, applies stringent filters to amplicon reads, and outputs regions that do not map to a reference sequence. It is aimed to build libraries of variants using kmers as a fast mapping tool, and to allow evaluation of the context in which variants are found. At this time, it has been set up and tested on fastq format sequences from Illumina, using 14mers. It is written in the language Go [or Golang] and takes about 6 seconds to run on a 2019 Mac Pro (3.5 GHz 8 core Intel Xeon W) for the test case of 5 Illumina sequence read files with filters already constructed. 
# Input and Output File Structure
AnCovMulti imports a number of standard packages as well as the custom globals and seqmerAC packages, which should be included in the release. 
The _**globals**_ package uses **control**, **mode**, and **factory** files to introduce environmental parameters, and these files must be included in the top-level folder. Parameters settings are in the form “_name = value_”, and if a required parameter is not in the mode or control file the program will barf. The program will also barf if the value cannot be converted to the expected type, but other than that there is no parameter checking. Parameters used in the program are _their last set value_, and _lines are not read_ past a _hashtag_. 
* The **control** file is the last parameter file read, and the parameter choices here will override any previous choices. Ideally it will be used to control most parameter changes among runs. 
* The **mode** file is meant to set the type of analysis being done, and the structure of the files. A lot of parameters are set in it. They are organized with headers. 
* The **factory** file is meant to set what the parameters, their type, and any bounds should be. It is meant to not be messed with and is not really implemented here. 
Using the naming conventions in the supplied mode file, 
* **sequences** are in fastq format with .fq ending, and are found in the sequences folder. 
* **Reference files** are in the references folder and include 
    -   the _**WuhanHu1**_ sequence in fasta format;

    -   _**left_primers.txt**_ and _**right_primers.txt**_ (tab-delimited start and
        stop in WuhanHu1 reference, no header);

    -   the _**WuhanWu1_kpos.xls**_ kmer, count, positions (tab-delimited,
        with header);

    -   _**SequenceDayTranslation.txt**_ in tab-delimited format;

    -   _**filterseqs_2020-10-14_CO.xls**_ list of filter names on different
        lines, not used but required;

    -   and _**kcounts14_WuhanHu1_14Oct2020.xls**_ which is not used and is
        the same as the kpos file without the positions.

-   **Query sequence filter** files, when produced, are placed into the
    qfilter folder. They are later read from the same folder.

-   All **other output** files are placed in the ACmulti_output folder,
    including

    -   _**Kpos.xls**_, the positions of the kmers in the reference sequence

        -   This should be the same as the kpos file in the refefences
            folder

    -   _**Primers_out.xls**_, first position, end position, direction (left
        or right) in tab-delimited format with headers.

    -   _**qnotk**_ contains qnotk files for the filtered query sequences,
        meaning they are kmers found in the query sequences but not in
        the reference sequence

        -   labeled left or right with klen and seqfile name, .xls at
            end, in tab-delimited format with headers.

        -   Kmer sequence, count in filtered query sequences, first kmer
            location in WuhanHu1 ref, read count (from pcounts file),
            calculated positions in Wuhan reference sequence based on
            primer location

    -   _**kandq**_ contains kandq files for the filtered query sequences,
        meaning they are kmers found in both the query sequences and the
        reference sequence

        -   labeled left or right with klen and seqfile name, .xls at
            end, in tab-delimited format with headers.

        -   Kmer sequence, count in filtered query sequences, first kmer
            location in WuhanHu1 ref, read count (from pcounts file),
            calculated positions in Wuhan reference sequence based on
            primer location

        -   Only printed if *bigoutput* parameter is true

    -   _**Kcounts**_ contains kcount files for the filtered query sequences,
        labeled left or right with klen and seqfile name, .xls at end,
        in tab-delimited format with headers.

        -   Kmer sequence, count in filtered query sequences, first kmer
            location in WuhanHu1 ref, read count (from pcounts file),
            calculated positions in Wuhan reference sequence based on
            primer location

        -   Contains all primers found in either qnotk or kandq

        -   Only printed if *bigoutput* parameter is true

    -   _**Pcounts**_ contain primer location and count found in filtered
        query sequence file

        -   This count is used to identify read counts

The provided folder structure contains the output expected with the
sequence, references, and control file parameter settings. (bigoutput is
false; qnotkmode is false but was run previously to produce the
qfilters)


# Control Parameters
## Variables likely to change in the control file

_qminprint = 1_
* don’t print qmers less than this to kcount file for the query and related kmer sets

_limitseqs = 50_
* just limit the number of seqs to look at for testing if you have a really big set of sequence files

_runorient = right_
* either right or left, to control the orientation of analysis

_qnotkmode = t_
* control run so don't repeat creation of qfilters
* t, T, true, True, TRUE accepted

_bigoutput= t_
* don't output big kcount and kandq files if we don't need them

## Variables not likely to change 

_klen = 14_
* klen needs to match the reference inputs or things will be a mess

_left = left_

_right = right_
* left and right need to match the reference input or it will be a mess

_printNs = false_
* do you really ever want to print kmers with Ns I them? Probably not. 

## Query filter boundary parameters  
_minseq = 247_	
* Illumina length minimum for 250 bp runs

_minmatch = 50_	
* the number of kmer matches in a read (matchcount) cannot be less than this

_maxmisses = 42_
* the number of misses is possible - matchcount + 1
* the firstmatch and lastmatch are the order of reference kmers found in the query sequence
* first and second are whichever of these are first and second in order in the Wuhan reference sequence

_Wudiffmax = 254_
* the difference between the second and first + klen cannot be more than this

_maxfirsthit = 2_
* the read position of the Wuhan first kmer cannot be greater than this

_minlasthit = 219_
* the read position of the Wuhan first kmer cannot be greater than this
* posdiff is the lastposition - sequence length + klen -1

_minposdiff = -16_
* posdiff can't be less than this

if left, the first kmer match must be a primer start

if right, the second kmer match + klen must be a primer end boundary


-   The remaining variables are less important and all variables are
    explained in the mode file

# Code
Three folders of [source code](https://github.com/PollockLaboratory/AnCovMulti/tree/main/src) are provided, the **main** code
* **[AnCovMulti](https://github.com/PollockLaboratory/AnCovMulti/tree/main/src/AnCovMulti)**

along with two support **packages**
1. **[seqmerAC](https://github.com/PollockLaboratory/AnCovMulti/tree/main/src/seqmerAC)** 
1. **[globals](https://github.com/PollockLaboratory/AnCovMulti/tree/main/src/globals)**

These contain, respectively, the files AnCovMulti.go, seqmerAC.go, and globals.go. 

The AnCovMulti folder also contains the [binary](https://github.com/PollockLaboratory/AnCovMulti/blob/main/src/AnCovMulti/AnCovMulti) AnCovMulti file, compiled on _macOS version 10.15.5_ and using _go1.14.4_ downloaded June 1 2020. 

Information about the Go programming language can be found at [Go](http://golang.org/). 

# Future
The near-term aim is to automate the variant analysis, which will output a list of variants of interest along with their position in the reference, what they changed from and into in nucleotides and amino acids (if in protein, also listed). Further, the filenames will be translated into human readable identifiers, such as day of sampling and source, and the frequencies will be provided across the files. There will also be a similar file with more details, such as left and right primer source, R1 and R2 for Illumina files, the mean count of the variant, or min and max, or best estimate, the number of reads of that read type, and the number of reference reads of that read type. Finally, we will address the Minion reads and figure out reasonable parameter settings for them.
