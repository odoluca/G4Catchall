# G4Catchall
G4Catchall is a python package designed to scan given DNA/RNA sequences for G-quadruplexes with or without atypical features

Please cite: 
Doluca, O. (2019). G4Catchall: A G-quadruplex prediction approach considering atypical features. Journal Of Theoretical Biology, 463, 92-98. doi: 10.1016/j.jtbi.2018.12.007

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
        echo ^>mychr > mychr.fa 
        echo TTGGGTTGGGACTGGGTACGGGAATAAATAGGTTAGGAATGGATAGGATCCCTTCCCTTCCCTTCCCTTGGCGCGGCCGGCGG >> mychr.fa
        

        python G4Catchall.py -f mychr.fa --G3L 1..3 
        mychr   2       22      20      +       GGGTTGGGACTGGGTACGGG    GGGTTGGGACTGGGTACGGG    1.7
        mychr   49      67      18      -       CCCTTCCCTTCCCTTCCC      GGGAAGGGAAGGGAAGGG      -2.0
        
        ## 2 Guanine-tetrad G-quadruplexes can be included using --G2L
        python G4Catchall.py -f mychr.fa --G3L 1..3 --G2L 1..3
        mychr   2       22      20      +       GGGTTGGGACTGGGTACGGG    GGGTTGGGACTGGGTACGGG    1.7
        mychr   30      47      17      +       GGTTAGGAATGGATAGG       GGTTAGGAATGGATAGG       0.9411764705882353
        mychr   49      67      18      -       CCCTTCCCTTCCCTTCCC      GGGAAGGGAAGGGAAGGG      -2.0

        ## Score threshold can be changed using --G4Threshold
        python G4Catchall.py -f mychr.fa --G3L 1..3 --G2L 1..3 --G4HThreshold 0.4
        mychr   2       22      20      +       GGGTTGGGACTGGGTACGGG    GGGTTGGGACTGGGTACGGG    1.7
        mychr   30      47      17      +       GGTTAGGAATGGATAGG       GGTTAGGAATGGATAGG       0.9411764705882353
        mychr   49      67      18      -       CCCTTCCCTTCCCTTCCC      GGGAAGGGAAGGGAAGGG      -2.0
        mychr   69      83      14      +       GGCGCGGCCGGCGG  GGCGCGGCCGGCGG  0.7142857142857143
                    
        ##When no fasta file is indicated, it only constructs the regex from given parameters and prints.
        python G4Catchall.py --G3L 1..3 -I 0
            ([Gg]{3,})  (\w{1,3})  ([Gg]{3,}) (\w{1,3}) ([Gg]{3,}) (\w{1,3}) ([Gg]{3,})
            
    DOWNLOAD
        G4Catchall.py is hosted at http://github.com/odoluca/G4Catchall
    
    PLEASE CITE
        Doluca, O. (2019). G4Catchall: A G-quadruplex prediction approach considering atypical features. Journal Of Theoretical Biology, 463, 92-98. doi: 10.1016/j.jtbi.2018.12.007

optional arguments:
  -h, --help            show this help message and exit
  --fasta FASTA, -f FASTA
                        Input file in fasta format containing one or more sequences can be used.
                        Please note that, if not used, only the regular expression constructed using given
                        arguments will be printed. 
                        
  --min_Gtract_for_extreme_loop MIN_GTRACT_FOR_EXTREME_LOOP, -E MIN_GTRACT_FOR_EXTREME_LOOP
                        Defines the minimum G-tract length for permission of an extreme 
                        loop. Works only with --extreme_loop. Can be set to 2 or 3. Default=3
                        
  --extreme_loop [EXTREME_LOOP], --XL [EXTREME_LOOP]
                        Allows search for an extreme loop. If precedes a secondary argument,
                        such as "1..20" also defines the limits of the loop. For default values do 
                        not use a second argument. Default="1..30"
                        
  --G2GQs_allowed       Allows G-quadruplexes with G-tracts of two guanines. Not necessary
                        with --G2GQ_loop command.
                        
  --G2GQ_loop [G2GQ_LOOP], --G2L [G2GQ_LOOP]
                        Allows G-quadruplexes with G-tracts of two guanines and defines 
                        limits of loops for such G-quadruplexes if precedes a secondary argument,
                        such as "1..7". Do not use a secondary argument for default loop limits.
                        Default="1..2" 
                        
  --G3GQ_loop G3GQ_LOOP, --G3L G3GQ_LOOP
                        Defines limits of loops for typical G-quadruplexes if precedes a 
                        secondary argument, such as "1..7". Do not use for default loop limits.
                        Default="1..8" 
                        
  --max_imperfect_Gtracts MAX_IMPERFECT_GTRACTS, -I MAX_IMPERFECT_GTRACTS
                        Defines the number of atypical or "imperfect" G-tracts allowed for
                        G-quadruplexes with G-tracts of at least 3 guanines. It can be set to 0,1 
                        or 2. Default=1
                        
  --bulge_only, -B      Defines the nature of the imperfect G-tracts allowed. If used only
                        bulged G-tracts are allowed. Otherwise, mismatches are also allowed.
                        
  --max_GQ_length MAX_GQ_LENGTH, --max MAX_GQ_LENGTH
                        Maximum allowed G-quadruplex length for a single discovery. This 
                        should be used with caution. If not used together with --dont_merge_overlapping
                        discovered sequences may be longer than the given value. This parameter is 
                        essentially designed for limiting cumulative negative impact of long loops.
                        
  --no_reverse, -R      By default the program searches both strands by reversing the regex. 
                        If used only + strand is searched for matches.
                        
  --dont_merge_overlapping
                        Putative G-quadruplex-forming sequences may be found overlapping 
                        on the same strand. By default the program merges these sequences. If 
                        used, these are not overlapped. Using may result in huge number of 
                        matches and cause memory issues.
                        
  --include_flanks, --F
                        By default the program extracts only matching sequences. 
                        If used flanking nucleotides are also included in the search. 
                        Please note if used G-quadruplex-forming sequences at the beginning or 
                        ending of the sequences may be missed. Consider adding "N" to the edges 
                        of the sequence if G-quadruplex forming sequences are expected to be 
                        found at the very edge of the target sequence.
                        
  --G4H                 By default the program extracts only matching sequences. If 
                        used, discovered sequences are evaluated based on G4Hunter algorithm.
                        Low G4Hunter scores can be eliminated using --G4HThreshold argument. 
                        
  --G4HThreshold G4HTHRESHOLD
                        Removes G-quadruplex predictions with lower scores than the preceding 
                        threshold value. If used, --G4H usage is not necessary.
