"""
DEKO_tool
"""
from crispeta_util import *

parser = argparse.ArgumentParser(description = "CRISPETa_pc is a variation of CRISPETa focused on the design pf CRISPR paired sgRNAs on desired annotated regions of a genome."+
                                               "Main diference with CRISPETa is that sgRNAs are searched only in exons of target genes."+
                                               "CRISPETa_pc uses as input a list of IDs and an annotation file (gtf format).")

parser.add_argument('-i', '--input',
                    required = True,
                    dest = 'infile',
                    action = 'store',
                    help = 'Path to input file with target names')

parser.add_argument('-o', '--output',
                    required = False,
                    dest = 'outfile',
                    action = 'store',
                    default = './',
                    help = 'Path and prefix for output files')

parser.add_argument('-gtf', '--gtf',
                    dest = 'gtf',
                    action = 'store',
                    required = True,
                    help = 'Genome annotation (GTF) file')

parser.add_argument('-appris', '--appris',
                    dest = 'appris',
                    action = 'store',
                    required = True,
                    help = 'Principal isoforms')

parser.add_argument('-g', '--genome',
                    dest = 'genome',
                    action = 'store',
                    required = True,
                    help = 'Genome sequence (fasta) file')

parser.add_argument('-c', '--construct_method',
                    dest = 'construct_method',
                    action = 'store',
                    default = None,
                    type = str,
                    choices=[None,'DECKO'],
                    help = 'Method applied when making sgRNA pairs and oligo construction')

parser.add_argument('-sm', '--scoring_method',
                    dest = 's_method',
                    action = 'store',
                    default = 'Doench',
                    choices=['Doench', 'Han'],
                    help = 'Scoring method')

parser.add_argument('-si', '--min_iscore',
                    dest = 'individual_score',
                    action = 'store',
                    default = None,
                    type= float,
                    help = 'Minimum individual sgRNA score allowed')

parser.add_argument('-sp', '--min_pscore',
                    dest = 'paired_score',
                    action = 'store',
                    default = None,
                    type = float,
                    help = 'Minimum sgRNA paired score allowed')

parser.add_argument('-t', '--off-targets',
                    dest = 'off_targets',
                    action = 'store',
                    default = '0,0,0,x,x',
                    help = 'Maximum number of off-targets allowed per mismatch')

parser.add_argument('-r', '--rank',
                    dest = 'rank',
                    action = 'store',
                    default = 'appearance',
                    choices=['appearance', 'dist', 'score'],
                    help = 'Ranking feature for pairs [score]')

parser.add_argument('-mind', '--min_dist',
                    dest = 'min_dist',
                    action = 'store',
                    default = 30,
                    help = 'Minimun number of bp between protospacers allowed [genomic positions].')

parser.add_argument('-maxd', '--max_dist',
                    dest = 'max_dist',
                    action = 'store',
                    default = 3500,
                    help = 'Maximum number of bp between protospacers allowed [genomic positions]. Value 0 dissable this option')

parser.add_argument('-n', '--n',
                    dest = 'n',
                    action = 'store',
                    default = 10,
                    help = 'Number of pairs returned')

parser.add_argument('-v', '--variety',
                    dest = 'variety',
                    action = 'store',
                    default = 0.2,
                    help = 'Max number of times a single sgRNA would appear on results')


start_time = datetime.datetime.now()

#Storing user options
options = parser.parse_args()
if options.individual_score == None:
     options.individual_score = 0.2 if options.s_method == 'Doench' else 0
if options.paired_score == None:
     options.paired_score = 0.4 if options.s_method == 'Doench' else 0

#Creating a temporary directory for temporary files and opening log file
#tmp_dir = tempfile.mkdtemp(prefix='sgRNA')

#extract ID's from input file
ids = get_IDs(options.infile)

#extract CCDS ID's
transcripts = get_transcript_ids(options.appris, ids)

#extract exons
regions = get_regions(options.gtf, transcripts)

#extract sequences
regions = get_sequences(regions, options.genome)

#extract sgRNAs
model = create_model()
sgrnas = get_sgRNAs(regions, options.individual_score, options.paired_score, options.off_targets, model, options.s_method)

#make pairs
pairs = get_pairs(sgrnas, options.min_dist, options.rank, options.max_dist)

#results
if options.outfile == "./":
     name = './CRISPETa_results'
else:
     name = options.outfile
file_browser = open(name+'.bed', 'w')
file_pairs = open(name+'.txt', 'w')

print_browser(file_browser, pairs, options.n, options.variety)

output_pairs = print_pairs(file_pairs, pairs, options.n, options.variety, options.construct_method)
print_log(ids, output_pairs, start_time, options)