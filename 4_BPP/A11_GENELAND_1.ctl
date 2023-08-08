         seed =  -1

       seqfile = genetic_matrix.txt
      Imapfile = Imap_GENELAND.txt
       outfile = out_geneland.txt
      mcmcfile = mcmc1.txt

* speciesdelimitation = 0 * fixed species tree
* speciesdelimitation = 1 0 2   * species delimitation rjMCMC algorithm0 and finetune(e)
 speciesdelimitation = 1 1 2 0.5  * species delimitation rjMCMC algorithm1 finetune (a m)
          speciestree = 1 * species-tree by NNI

  speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted

  species&tree = 4  A  B  C  D 
                   78  44  8  4  
                  (B,(D,(A,C)));
        phase =   0  0  0  0 
                  
       usedata = 1  * 0: no data (prior); 1:seq like
         nloci = 4  * number of data sets in seqfile

     cleandata = 1    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = 3 0.02  # invgamma(a, b) for theta
      tauprior = 3 0.001    # invgamma(a, b) for root tau & Dirichlet(a) for other tau's

     heredity = 2 heredity.txt
     locusrate = 1 0 0 5 iid # (1: estimate locus rates mui, mubar=1 fixed)

      finetune = 1: 5 0.001 0.001  0.001 0.3 0.33   # finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars, Genetrees
        burnin = 8000
      sampfreq = 2
       nsample = 100000
