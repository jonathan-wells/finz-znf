  seqfile = /Users/jonwells/Projects/feschottelab/finz-znf/data/selection-analysis/rerio_subgroup.fa         * sequence data filename
 treefile = /Users/jonwells/Projects/feschottelab/finz-znf/data/selection-analysis/rerio_subgroup.fa.treefile
  outfile = rerio_subgroup.out   * main result file name

    noisy = 3      * 0,1,2,3,9: how much rubbish on the screen
  verbose = 1      * 1:detailed output
  runmode = 0      * 0:user-defined tree

  seqtype = 1      * 1:codons
CodonFreq = 1      * 0:equal, 1:F1X4, 2:F3X4, 3:F61
    model = 0      * 0: Bran
  NSsites = 1 2    * 
    icode = 0      * 0:universal code

fix_kappa = 1      * 1:kappa fixed, 0:kappa to be estimated
    kappa = 2.379  * initial or fixed kappa

fix_omega = 0      * 1:omega fixed, 0:omega to be estimated 
    omega = 0.751      *  
