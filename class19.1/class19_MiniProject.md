# class19_MiniProject
Leah Johnson, PID: A17394690

From a BLASTp search against human refseq: Official symbol: MYC Official
full name: MYC proto-oncogene protein isoform 1, bHLH transcription
factor \[Homo Sapiens\]

``` r
library(bio3d)

fast_seq <- read.fasta("A17394690_mutant_seq.fa")
fast_seq
```

                   1        .         .         .         .         .         60 
    wt_healthy     MDFFRVVENQQPPATMPLNVSFTNRNYDLDYDSVQPYFYCDEEENFYQQQQQSELQPPAP
    mutant_tumor   MDFFRVVENQQPPATMPLNVSFTNRNYDLDYDSVQPYFYCDEEENFYQQQQQSELQPPAP
                   ************************************************************ 
                   1        .         .         .         .         .         60 

                  61        .         .         .         .         .         120 
    wt_healthy     SEDIWKKFELLPTPPLSPSRRSGLCSPSYVAVTPFSLRGDNDGGGGSFSTADQLEMVTEL
    mutant_tumor   SEDIWKKFELLPTPPLSPSRRSGLVSPSYVAVTPFSLRGDNDGGGGSFSTADQLEMVTEL
                   ************************ *********************************** 
                  61        .         .         .         .         .         120 

                 121        .         .         .         .         .         180 
    wt_healthy     LGGDMVNQSFICDPDDETFIKNIIIQDCMWSGFSAAAKLVSEKLASYQAARKDSGSPNPA
    mutant_tumor   LGGDMVNQSFICDPDDETFIKNIIIQDCMWSGFSEAAKLVSEKLASYQRARKDSGSPNPA
                   ********************************** ************* *********** 
                 121        .         .         .         .         .         180 

                 181        .         .         .         .         .         240 
    wt_healthy     RGHSVCSTSSLYLQDLSAAASECIDPSVVFPYPLNDSSSPKSCASQDSSAFSPSSDSLLS
    mutant_tumor   RGHSVCSTSSLYLQDLSAAASECIDPSVVFPYPLNDSSSPKSCASQDSSAFSPSSDSLLS
                   ************************************************************ 
                 181        .         .         .         .         .         240 

                 241        .         .         .         .         .         300 
    wt_healthy     STESSPQGSPEPLVLHEETPPTTSSDSEEEQEDEEEIDVVSVEKRQAPGKRSESGSPSAG
    mutant_tumor   STESSPQGSPEPLVLHEETPPTTSSDSEEEQEDEYEIDVVSVEKRQAPGKRSESGSPSAG
                   ********************************** ************************* 
                 241        .         .         .         .         .         300 

                 301        .         .         .         .         .         360 
    wt_healthy     GHSKPPHSPLVLKRCHVSTHQHNYAAPPSTRKDYPAAKRVKLDSVRVLRQISNNRKCTSP
    mutant_tumor   GHSKPPHSPLVLKRCHVSTHQHNYAAPPSTRKDYPAAKRVKLDSVRVLRQISNNRKCTSP
                   ************************************************************ 
                 301        .         .         .         .         .         360 

                 361        .         .         .         .         .         420 
    wt_healthy     RSSDTEENVKRRTHNVLERQRRNELKRSFFALRDQIPELENNEKAPKVVILKKATAYILS
    mutant_tumor   RSSDTEENVKRRTHNVLERQRRNELKRSFFALRDQIPELENNEKAPKVVILKKATAYILS
                   ************************************************************ 
                 361        .         .         .         .         .         420 

                 421        .         .         .   454 
    wt_healthy     VQAEEQKLISEEDLLRKRREQLKHKLEQLRNSCA
    mutant_tumor   VQAEEQKLISEEDLLRKRREQLKHKLEQLRNSCA
                   ********************************** 
                 421        .         .         .   454 

    Call:
      read.fasta(file = "A17394690_mutant_seq.fa")

    Class:
      fasta

    Alignment dimensions:
      2 sequence rows; 454 position columns (454 non-gap, 0 gap) 

    + attr: id, ali, call

We can score residue conservation then find the non 1.0 scoring
positions. These will be the mutation positions.

``` r
conserv(fast_seq)
```

      [1]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
     [16]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
     [31]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
     [46]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
     [61]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
     [76]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0 -0.2  1.0  1.0  1.0  1.0  1.0
     [91]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [106]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [121]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [136]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [151]  1.0  1.0  1.0  1.0 -0.1  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [166]  1.0  1.0  1.0 -0.1  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [181]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [196]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [211]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [226]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [241]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [256]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [271]  1.0  1.0  1.0  1.0 -0.4  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [286]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [301]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [316]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [331]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [346]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [361]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [376]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [391]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [406]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [421]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [436]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [451]  1.0  1.0  1.0  1.0

``` r
mutation.sites <- which(conserv(fast_seq)<1)
```

``` r
paste(fast_seq$ali[1,mutation.sites], 
      mutation.sites, 
      fast_seq$ali[2,mutation.sites])
```

    [1] "C 85 V"  "A 155 E" "A 169 R" "E 275 Y"
