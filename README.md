# G4Catchall
G4Catchall is a python package designed to scan given DNA/RNA sequences for G-quadruplexes with or without atypical features

PUBLICATION IN REVIEW...

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
        
