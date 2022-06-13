table bigPrimerBed
"BigBed schema for Primer JuJu primer pairs"
    (
    string chrom;      "Chromosome"
    uint   chromStart; "Start position in chromosome"
    uint   chromEnd;   "End position in chromosome"
    string name;       "Name of item"
    uint   score;      "Score from 0-1000"
    char[1] strand;    "+ or -"
    uint thickStart;   "Start of where display should be thick"
    uint thickEnd;     "End of where display should be thick"
    uint reserved;     "itemRgb"
    int blockCount;    "Number of blocks"
    int[blockCount] blockSizes; "Comma separated list of block sizes"
    int[blockCount] chromStarts; "Start positions relative to chromStart"
    int transcritome_off_target_cnt  "transcriptome off-target hits"
    int genome_off_target_cnt  "genomic off-target hits"
    int    PRIMER_PAIR_PRODUCT_SIZE;      "PAIR_PRODUCT_SIZE"
    float  PRIMER_PAIR_COMPL_ANY_TH;      "PAIR_COMPL_ANY_TH"
    float  PRIMER_PAIR_COMPL_END_TH;      "PAIR_COMPL_END_TH"
    float  PRIMER_PAIR_PENALTY;           "PAIR_PENALTY"
    string PRIMER_LEFT;                   "LEFT"
    string PRIMER_LEFT_SEQUENCE;          "LEFT_SEQUENCE"
    float  PRIMER_LEFT_END_STABILITY;     "LEFT_END_STABILITY"
    float  PRIMER_LEFT_GC_PERCENT;        "LEFT_GC_PERCENT"
    float  PRIMER_LEFT_HAIRPIN_TH;        "LEFT_HAIRPIN_TH"
    float  PRIMER_LEFT_PENALTY;           "LEFT_PENALTY"
    float  PRIMER_LEFT_SELF_ANY_TH;       "LEFT_SELF_ANY_TH"
    float  PRIMER_LEFT_SELF_END_TH;       "LEFT_SELF_END_TH"
    float  PRIMER_LEFT_TM;                "LEFT_TM"
    string PRIMER_RIGHT;                  "RIGHT"
    string PRIMER_RIGHT_SEQUENCE;         "RIGHT_SEQUENCE"
    float  PRIMER_RIGHT_END_STABILITY;    "RIGHT_END_STABILITY"
    float  PRIMER_RIGHT_GC_PERCENT;       "RIGHT_GC_PERCENT"
    float  PRIMER_RIGHT_HAIRPIN_TH;       "RIGHT_HAIRPIN_TH"
    float  PRIMER_RIGHT_PENALTY;          "RIGHT_PENALTY"
    float  PRIMER_RIGHT_SELF_ANY_TH;      "RIGHT_SELF_ANY_TH"
    float  PRIMER_RIGHT_SELF_END_TH;      "RIGHT_SELF_END_TH"
    float  PRIMER_RIGHT_TM;               "RIGHT_TM"
    )
