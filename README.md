This perl script takes a list of miRNA sequences (in fasta format) and searches for near-perfect complements in the transcript sequence (also in fasta format) using program scan_for_matches, evaluates the miRNA:target complement patterns, calculate minimum free energie ratio (MFE of target:miRNA / MFE of miRNA perfect match) using program RNAhybrid, and outputs the prediction results in an easy-to-read html web table format.

External programs:
scan_for_matches: https://github.com/hallundbaek/scan_for_matches.git

RNAhybrid: https://bibiserv.cebitec.uni-bielefeld.de/rnahybrid

Please cite: Sun YH, Lu S, Shi R, Chiang VL. Computational prediction of plant miRNA targets. Methods Mol Biol. 2011;744:175-86. doi: 10.1007/978-1-61779-123-9_12. PMID: 21533693.
