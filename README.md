# CRISPETa_PC


### Example:
##### Input file [-i]: Gene name according to GENCODE, one name per line
$ cat input.txt
  SLC2A1
  HK1
  PLIN2
  AKT1
  PTEN

##### Running command:
$ python CRISPETa_coding.py -i input.txt \
			    -gtf human.gtf \
			    -appris appris.txt \
			    -g human.fa

