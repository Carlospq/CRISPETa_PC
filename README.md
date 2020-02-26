# CRISPETa_PC
Adaptation of CRISPETa for paired sgRNA designs in protein coding ORFs\
CRISPETa_PC needs a database (sqlite) in order to calculate the number of Off-Targets. Check reference for more information


<i>Pulido-Quetglas, C., Aparicio-Prat, E., Arnan, C., Polidori, T., Hermoso, T., Palumbo, E., Ponomarenko, J., Guigo, R., Johnson, A. K. (2017). Scalable Design of Paired CRISPR Guide RNAs for Genomic Deletion. PLOS Computational Biology, 13(3), e1005341. https://doi.org/10.1371/journal.pcbi.1005341


### Example:
##### Input file [-i]: Gene name according to GENCODE, one name per line
$ cat input.txt\
&nbsp;&nbsp;SLC2A1\
&nbsp;&nbsp;HK1\
&nbsp;&nbsp;PLIN2\
&nbsp;&nbsp;AKT1\
&nbsp;&nbsp;PTEN

##### Running command:
$ python CRISPETa_coding.py -i input.txt\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-gtf human.gtf\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-appris appris.txt\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-g human.fa

