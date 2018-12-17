from __future__ import division,print_function,absolute_import
import argparse, sys, re, string, regex



" --------------------------[ Parse arguments ]---------------------- "

parser=argparse.ArgumentParser(
    description="""
    DESCRIPTION
    
        Searches for matches to a G-quadruplex-fitting regex in a fasta file, 
        filters through G4Hunter-like secondary scoring scheme and return a bed file with 
        coordinates of the match, matched sequence, G-quadruplex forming sequence and the score.
        
        Output bed file has the following columns:
        1. description of the fasta sequence (e.g. NC_00024.11 Y chromosome)
        2. start of the match
        3. end of the match
        4. size of the match
        5. strand of the match (e.g. +)
        6. positive strand sequence of the match (e.g. CCCTTCCCTTTCCCTCCC)
        7. matched G-quadruplex-forming sequence (e.g. GGGAGGGAAAGGGAAGGG)
        8. score of the matched G-quadruplex-forming sequence based on selected scoring scheme
        
    EXAMPLE
        ##Test data:
        echo '>mychr' > /tmp/mychr.fa 
        echo 'TTGGGTTGGGACTGGGTACGGGAATAAATAGGTTAGGAATGGATAGGAT' >> /tmp/mychr.fa
        
        G4Catchall.py -f /tmp/mychr.fa --G3L 1..3 --G2L 1..3
            mychr   2   22  20  +   GGGTTGGGACTGGGTACGGG    GGGTTGGGACTGGGTACGGG
            mychr   30  47  17  +   GGTTAGGAATGGATAGG   GGTTAGGAATGGATAGG
        
        ## G4Hunter scores can be calculated and included at the end of the line.
        G4Catchall.py -f /tmp/mychr.fa --G3L 1..3 --G4H
            mychr   2   22  20  +   GGGTTGGGACTGGGTACGGG    GGGTTGGGACTGGGTACGGG    2.11
                    
        ##When no fasta file is indicated, it only constructs the regex from given parameters and prints.
        G4Catchall.py -I 0
            ([Gg]{3,})  (\w{1,8})  ([Gg]{3,}) (\w{1,8}) ([Gg]{3,}) (\w{1,8}) ([Gg]{3,})
        
    DOWNLOAD
        G4Catchall.py is hosted at http://github.com/odoluca/G4Catchall

""",
    formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--fasta', '-f',
                    type=str,
help='''Input file in fasta format containing one or more sequences can be used.
Please note that, if not used, only the regular expression constructed using given
arguments will be printed. 

'''
                    )

parser.add_argument('--min_Gtract_for_extreme_loop','-E',
type=int,
help="""Defines the minimum G-tract length for permission of an extreme 
loop. Works only with --extreme_loop. Can be set to 2 or 3. Default=3

""", default=3)

parser.add_argument('--extreme_loop','--XL',
type=str,nargs='?',const='1..30',
help="""Allows search for an extreme loop. If precedes a secondary argument,
such as "1..20" also defines the limits of the loop. For default values do 
not use a second argument. Default="1..30"

""", default=False)

parser.add_argument('--G2GQs_allowed',
action='store_true',
help="""Allows G-quadruplexes with G-tracts of two guanines. Not necessary
with --G2GQ_loop command.

""", default=False)

parser.add_argument('--G2GQ_loop','--G2L',
type=str,nargs='?',const='1..2',
help="""Allows G-quadruplexes with G-tracts of two guanines and defines 
limits of loops for such G-quadruplexes if precedes a secondary argument,
such as "1..7". Do not use a secondary argument for default loop limits.
Default="1..2" 

""", default=False)

parser.add_argument('--G3GQ_loop','--G3L',
type=str,
help="""Defines limits of loops for typical G-quadruplexes if precedes a 
secondary argument, such as "1..7". Do not use for default loop limits.
Default="1..8" 

""", default='1..8')

parser.add_argument('--max_imperfect_Gtracts','-I',
type=int,
help="""Defines the number of atypical or "imperfect" G-tracts allowed for
G-quadruplexes with G-tracts of at least 3 guanines. It can be set to 0,1 
or 2. Default=1

""", default=1)

parser.add_argument('--bulge_only','-B',
action='store_true',
help="""Defines the nature of the imperfect G-tracts allowed. If used only
bulged G-tracts are allowed. Otherwise, mismatches are also allowed.

""", default=False)

parser.add_argument('--max_GQ_length','--max',
                    type=int,
help='''Maximum allowed G-quadruplex length for a single discovery. This 
should be used with caution. If not used together with --dont_merge_overlapping
discovered sequences may be longer than the given value. This parameter is 
essentially designed for limiting cumulative negative impact of long loops.

''', default=False
                    )


parser.add_argument('--no_reverse','-R',
action='store_true',
help="""By default the program searches both strands by reversing the regex. 
If used only + strand is searched for matches.

""")

parser.add_argument('--dont_merge_overlapping',
action='store_true',
help="""Putative G-quadruplex-forming sequences may be found overlapping 
on the same strand. By default the program merges these sequences. If 
used, these are not overlapped. Using may result in huge number of 
matches and cause memory issues.

""")

parser.add_argument('--include_flanks','--F',
action='store_true',
help="""By default the program extracts only matching sequences. 
If used flanking nucleotides are also included in the search. 
Please note if used G-quadruplex-forming sequences at the beginning or 
ending of the sequences may be missed. Consider adding "N" to the edges 
of the sequence if G-quadruplex forming sequences are expected to be 
found at the very edge of the target sequence.

""")

parser.add_argument('--G4H',
action='store_true',
help="""By default the program extracts only matching sequences. If 
used, discovered sequences are evaluated based on G4Hunter algorithm.
Low G4Hunter scores can be eliminated using --G4HThreshold argument. 

""")

parser.add_argument('--G4HThreshold',
type=float,
help="""Removes G-quadruplex predictions with lower scores than the preceding 
threshold value. If used, --G4H usage is not necessary.

""")

args=parser.parse_args()

" --------------------------[ End of Parse arguments ]---------------------- "


" ------------------------------[  Reverse complement ]--------------------------------- "
def ReverseComplement(seq):
    seq1 = 'ATCGNTAGCNatcgntagcn'
    seq_dict = {seq1[i]: seq1[i + 5] for i in range(20) if i < 5 or 10 <= i < 15}
    return "".join([seq_dict[base] for base in reversed(seq)])

" ------------------------------[  End of Reverse complement  ]--------------------------------- "


""" ------------------------------[  Sorter ]--------------------------------- 
Code to sort list of lists
see http://www.saltycrane.com/blog/2007/12/how-to-sort-table-by-columns-in-python/
"""

def sort_table(table, cols):
    """ sort a table by multiple columns
        table: a list of lists (or tuple of tuples) where each inner list
               represents a row
        cols:  a list (or tuple) specifying the column numbers to sort by
               e.g. (1,0) would sort by column 1, then by column 0
    """
    for col in reversed(cols):
        table = sorted(table, key=operator.itemgetter(col))
    return (table)

" --------------------------[ End of Sorter ]---------------------- "


" ------------------------------[ Parameters  ]--------------------------------- "

" ------------------------------[ Parse Arguments ]--------------------------------- "

G4H_scoring=False
max_G4H_score=4


if args.G4HThreshold is not None: args.G4HunterScores=True

G2sAllowed =False
if args.G2GQs_allowed: G2sAllowed=True
shrtLoopMin,shrtLoopMax='1','2'
if args.G2GQ_loop:
    G2sAllowed = True
    shrtLoopMin,shrtLoopMax= args.G2GQ_loop.split("..")
typLoopMin,typLoopMax=args.G3GQ_loop.split("..")
extLoopMin, extLoopMax = '1','30'
ExtremeAllowed=False #args.extreme_loopAllowed
ExtremeAllowedForG2s = False
if args.extreme_loop:
    ExtremeAllowed=True
    extLoopMin,extLoopMax=args.extreme_loop.split("..")
if args.min_Gtract_for_extreme_loop==2:
    ExtremeAllowedForG2s=True
elif args.min_Gtract_for_extreme_loop==3:
    ExtremeAllowedForG2s=False
ImperfectTractsAllowed=args.max_imperfect_Gtracts
BulgedTractsOnly=args.bulge_only


if not G2sAllowed or not ExtremeAllowed:
    ExtremeAllowedForG2s = False

InclFlanks=args.include_flanks
if args.dont_merge_overlapping: MergeOverlapping=False
else: MergeOverlapping=True
NoReverse=args.no_reverse



" ------------------------------[  Define RegEx  ]--------------------------------- "


bulgeOnly="[Gg]{2,}[ATUCatuc][Gg]+|[Gg]+[ATUCatuc][Gg]{2,}"
mismatch="[Gg]{2,}|[Gg]+[ATUCatuc][Gg]+"
Dimp1='?P<imp1>'
Dimp2='?P<imp2>'
imp='('+mismatch+')'
if BulgedTractsOnly:
    imp='('+bulgeOnly+')'
Timp1='?(imp1)'
Timp2='?(imp2)'
shrt='\w{'+shrtLoopMin+','+shrtLoopMax+'}'
Tshrt='?(G2GQ)'
Dshrt='?P<G2GQ>[Gg]{2}'
typ='\w{'+typLoopMin+','+typLoopMax+'}'
ext='\w{'+extLoopMin+','+extLoopMax+'}'
Dext='?P<extLoop>'
Text='?(extLoop)'

# Construct Tract 1:
Tract1='[Gg]{3,}'
if ImperfectTractsAllowed>0: Tract1=Tract1+'|('+Dimp1+imp+')'
if G2sAllowed: Tract1=Tract1+'|('+Dshrt+')'

# Construct Loop 1:
Loop1=typ
if ExtremeAllowed:Loop1=Loop1+'|('+Dext+ext+')'
if G2sAllowed:
    shrtAdd = shrt
    if ExtremeAllowedForG2s:
        shrtAdd = shrtAdd + '|(' + Dext + ext+')'
    Loop1 = Tshrt + '(' + shrtAdd + ')|('+Loop1+')'

# Construct Tract 2:
Tract2='[Gg]{3,}'
if ImperfectTractsAllowed>1: Tract2=Tract2+'|('+Dimp2+imp+')'
if ImperfectTractsAllowed>0: Tract2=Timp1+'('+Tract2+')|([Gg]{3,}|('+Dimp1+imp+'))'
if G2sAllowed:
    Tract2=Tshrt+'[Gg]{2,}|('+Tract2+')'

# Construct Loop 2:
Loop2=typ
if ExtremeAllowed:
    Loop2=Text+Loop2+'|('+typ+'|('+Dext+ext+'))'
if G2sAllowed:
    shrtAdd=shrt
    if ExtremeAllowedForG2s:
        shrtAdd = Text+ shrtAdd + '|('+shrt+'|(' + Dext + ext+'))'
    Loop2=Tshrt+'('+shrtAdd+')|('+Loop2+')'

# Construct Tract 3:
Tract3='[Gg]{3,}'
if ImperfectTractsAllowed>1:
    Tract3=Timp2+Tract3+'|([Gg]{3,}|('+Dimp2+imp+'))'
if ImperfectTractsAllowed>0:
    Tract3=Timp1+'('+Tract3+')|([Gg]{3,}|('+Dimp1+imp+'))'
if G2sAllowed:
    Tract3=Tshrt+'[Gg]{2,}|('+Tract3+')'

# Combine all regions:
reg=r'('+Tract1+')  ('+Loop1+')  ('+ Tract2+') ('+Loop2+') ('+Tract3+') ('+Loop2+') ('+Tract3+')'
if InclFlanks: reg=r'\w'+reg+r'\w'


"""if no fasta is used, only return the regex and exit"""

if args.fasta is None:
    print (reg)
    exit()

""" Reverse forward match """
intab = 'actguACTGU'
outtab = 'tgacaTGACA'

if args.no_reverse is False:
    transtab = string.maketrans(intab, outtab)
    regrev = reg.translate(transtab)
else:
    regrev = ''



" ------------------------------[ End of Define RegEx ]--------------------------------- "


" ------------------------------[ Check Fasta ]--------------------------------- "

if args.fasta == '-':
    args.fasta = sys.stdin.readlines()
    if len(args.fasta) > 1:
        sys.exit('\nquadpareser.py: Only one input file at a time can be processed:\n--fasta/-f: %s\n' % (args.fasta))
    args.fasta = args.fasta[0].strip()
" ------------------------------[ End of Check Fasta ]--------------------------------- "


" ------------------------------[ Functions ]--------------------------------- "

""" Reverse complement function """
def ReverseComplement(seq):
    seq1 = 'ATCGNWSMKRYBDHVatcgnwsmkrybdhv'
    seq2 = 'TAGCNWSKMYRVHDBtagcnwskmyrvhdb'
    seq_dict = {seq1[i]: seq2[i] for i in range(len(seq1))}
    return "".join([seq_dict[base] for base in reversed(seq)])

"""                               G4Hunter
Modified code to score the discovered sequences based on G4Hunter algorithm. 
see Amina Bedrat, Laurent Lacroix, Jean-Louis Mergny; Re-evaluation of G-quadruplex propensity
 with G4Hunter, Nucleic Acids Research, Volume 44, Issue 4, 29 February 2016, Pages 1746-1759,
https://doi.org/10.1093/nar/gkw006
"""

def G4HScore(seq,minRepeat=2,penalizeGC=True):
    i=0
    baseScore=[]
    while i<len(seq):
        tractScore=[0]
        k=1
        GTract=False
        while seq[i]=="G":
            tractScore=[(min(k,4))] #derivation from original algorithm: tractScore=[min(k-1,16)]*k
            # region derivation from original algorithm: if prev is "C" apply bigger penalty. penalizes GCs
            if penalizeGC:
                try:
                    pass
                    if seq[i-k]=="C":baseScore[-1]=-2
                except:
                    pass
            # endregion
            k+=1
            i+=1
            GTract=True
            if i==len(seq): break
        if not GTract:
            while seq[i]=="C":
                tractScore=[max(-k,-4)] #derivation from original algorithm: tractScore=[max(-k,-16)]*k
                # region derivation from original algorithm: if prev is "G" apply bigger penalty. penalizes GCs
                if penalizeGC:
                    try:
                        pass
                        if seq[i - k] == "G": baseScore[-1] = 2
                    except:
                        pass
                # endregion
                k+=1
                i+=1
                GTract=True
                if i == len(seq): break
        baseScore=baseScore.__add__(tractScore)
        if not GTract: i += 1
        # print baseScore
    Score=0
    for value in baseScore:
        Score+=value
    return float(Score)/len(seq)

"""                               LIST SORTER
Code to sort list of lists
see http://www.saltycrane.com/blog/2007/12/how-to-sort-table-by-columns-in-python/
"""
import operator

def sort_table(table, cols):
    """ sort a table by multiple columns
        table: a list of lists (or tuple of tuples) where each inner list
               represents a row
        cols:  a list (or tuple) specifying the column numbers to sort by
               e.g. (1,0) would sort by column 1, then by column 0
    """
    for col in reversed(cols):
        table = sorted(table, key=operator.itemgetter(col))
    return (table)

" ------------------------------[ End of Functions ]--------------------------------- "



psq_re_f = regex.Regex(reg,regex.VERBOSE|regex.MULTILINE)
psq_re_r = regex.Regex(regrev,regex.VERBOSE|regex.MULTILINE)

ref_seq_fh = open(args.fasta)

ref_seq = []
line = (ref_seq_fh.readline()).strip()
chr = re.sub('^>', '', line)
line = (ref_seq_fh.readline()).strip()
gquad_list = []
while True:
    while line.startswith('>') is False:
        ref_seq.append(line)
        line = (ref_seq_fh.readline()).strip()
        if line == '':
            break
    ref_seq = ''.join(ref_seq)
    for m in psq_re_f.finditer(ref_seq,overlapped=True):
        if not args.max_GQ_length or len(m.group(0))<=args.max_GQ_length:
            if MergeOverlapping and (len(gquad_list)>0) and gquad_list[-1][4]=="+" and (m.start()<=gquad_list[-1][2]) and (chr==gquad_list[-1][0]):
                orj=gquad_list[-1]
                new_seq=orj[5]+m.group(0)[orj[2]-m.start():]
                G4Hscore=G4HScore(new_seq,2,True)
                if abs(G4Hscore)>=(args.G4HThreshold):
                    gquad_list[-1]=[chr, orj[1], m.end(), m.end()-orj[1], '+', new_seq,
                                new_seq,G4HScore(new_seq,2,True)]
            else:
                G4Hscore=G4HScore(m.group(0),2,True)
                if abs(G4Hscore)>=(args.G4HThreshold):
                    gquad_list.append([chr, m.start(), m.end(), len(m.group(0)), '+', m.group(0),
                                 m.group(0),G4HScore(m.group(0),2,True)])  # modification: added sequence again
    if args.no_reverse is False:
        if not args.max_GQ_length or len(m.group(0)) <= args.max_GQ_length:
            for m in psq_re_r.finditer(ref_seq,overlapped=True):
                if MergeOverlapping and (len(gquad_list) > 0) and gquad_list[-1][4]=="-" and (m.start() <= gquad_list[-1][2]) and (chr == gquad_list[-1][0]):
                    orj = gquad_list[-1]
                    new_seq = orj[5] + m.group(0)[orj[2] - m.start():]
                    G4Hscore = G4HScore(new_seq, 2, True)
                    if abs(G4Hscore) >= (args.G4HThreshold):
                        gquad_list[-1] = [chr, orj[1], m.end(), m.end() - orj[1], '-', new_seq,
                                      ReverseComplement(new_seq),G4HScore(new_seq,2,True)]
                else:
                    G4Hscore = G4HScore(m.group(0), 2, True)
                    if abs(G4Hscore) >= (args.G4HThreshold):
                        gquad_list.append([chr, m.start(), m.end(), len(m.group(0)), '-', m.group(0),
                                       ReverseComplement(m.group(0)), G4HScore(m.group(0),2,True)])  # modification: added reverse complement


    chr = re.sub('^>', '', line)
    ref_seq = []
    line = (ref_seq_fh.readline()).strip()
    if line == '':
        break

gquad_sorted = sort_table(gquad_list, (0, 1, 2, 3))



for line in gquad_sorted:
    line = '\t'.join([str(x) for x in line])
    print (line)
sys.exit()
