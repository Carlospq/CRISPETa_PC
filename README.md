# CRISPETa_PC


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

