# AnCovMulti
Ancestral Covid Variant Analysis

The current program v0.1 created Dec. 28, 2020

AnCovMulti is written in Golang to analyze variants in Covid amplicon sequence data. The first version was written December 28, 2020 by David Pollock with consultation from Richard Goldstein. It is currently restricted to Illumina reads. 
Please cite Kemp et al., 2020*.

Further information about the program and how to run it can be found in the [Tutorial](https://github.com/PollockLaboratory/AnCovMulti/wiki/Tutorial) or AnCovMulti_Readme3.docx

All code is to be found in [source code](https://github.com/PollockLaboratory/AnCovMulti/tree/main/src)

An [Example](https://github.com/PollockLaboratory/AnCovMulti/tree/main/Example) includes a set of reference files and a small sequence file. Note that the example [sequence file](NB13_CAMB-1B5CE7_L001_R1.fq.zip) is compressed because it is 50Mb when unzipped, but it runs quickly and you need at least a single fastq file from a sequencing run to apply the program sensibly. When uncompressed this should be deposited in the folder named "sequences".

The behaviors if you deviate from the prescribed file and folder formats are unknown.

#References

**SARS-CoV-2 evolution during treatment of chronic infection**. [pubmed](https://pubmed.ncbi.nlm.nih.gov/33398302/)
Kemp SA, Collier DA, Datir RP, Ferreira IATM, Gayed S, Jahun A, Hosmillo M, Rees-Spear C, Mlcochova P, Lumb IU, Roberts DJ, Chandra A, Temperton N; CITIID-NIHR BioResource COVID-19 Collaboration; COVID-19 Genomics UK (COG-UK) Consortium, Sharrocks K, Blane E, Modis Y, Leigh KE, Briggs JAG, van Gils MJ, Smith KGC, Bradley JR, Smith C, Doffinger R, Ceron-Gutierrez L, Barcenas-Morales G, Pollock DD, Goldstein RA, Smielewska A, Skittrall JP, Gouliouris T, Goodfellow IG, Gkrania-Klotsas E, Illingworth CJR, McCoy LE, Gupta RK. _**Nature**_. 2021 Apr;592(7853):277-282. doi: 10.1038/s41586-021-03291-y. Epub 2021 Feb 5. PMID: 33545711 

**Neutralising antibodies in Spike mediated SARS-CoV-2 adaptation**. [pubmed](https://pubmed.ncbi.nlm.nih.gov/33545711/)
Kemp SA, Collier DA, Datir R, Ferreira I, Gayed S, Jahun A, Hosmillo M, Rees-Spear C, Mlcochova P, Lumb IU, Roberts DJ, Chandra A, Temperton N; CITIID-NIHR BioResource COVID-19 Collaboration; COVID-19 Genomics UK (COG-UK) Consortium, Sharrocks K, Blane E, Briggs J, van GM, Smith K, Bradley JR, Smith C, Doffinger R, Ceron-Gutierrez L, Barcenas-Morales G, Pollock DD, Goldstein RA, Smielewska A, Skittrall JP, Gouliouris T, Goodfellow IG, Gkrania-Klotsas E, Illingworth C, McCoy LE, Gupta RK. _**medRxiv**_. 2020 Dec 29:2020.12.05.20241927. doi: 10.1101/2020.12.05.20241927. Preprint. PMID: 33398302 






