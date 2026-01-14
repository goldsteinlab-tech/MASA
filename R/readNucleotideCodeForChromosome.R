
#' Read or Generate Nucleotide Codes for a Chromosome
#'
#' Reads or computes nucleotide codes for a specified chromosome. If a precomputed code file exists, it loads the codes from disk; otherwise, it computes the codes from the genome sequence, optionally masking unmappable bases. Returns a list containing the nucleotide codes and a frequency table.
#'
#' @param chr Character. Chromosome name or identifier.
#' @param nmer Integer. Length of the nucleotide sequence to consider (default is 4).
#' @param genome Named list or object. Genome sequence containing chromosome data.
#' @param nuccodefileDir Character. Directory where nucleotide code files are stored (default is current directory).
#' @param mapdir Character. Directory where mappability files are stored (default is empty string).
#'
#' @return List with two elements:
#'   \item{code}{Integer vector of nucleotide codes for the chromosome.}
#'   \item{freqtable}{Data frame with columns \code{seq} (nucleotide sequence label) and \code{count} (frequency of each code).}
#' The result is also saved as an RDS file in \code{nuccodefileDir}.
#' @export
#'
#' @examples
#' # Example usage:
#' # genome <- list(chr1 = "ACGTACGTACGT")
#' # result <- readNucleotideCodeForChromosome("chr1", nmer = 4, genome = genome)
readNucleotideCodeForChromosome<-function(chr, nmer=4, genome=NA, nuccodefileDir='.', mapdir='') {

  unmappable = NA;

  np =nmer;
  #nuccodefile = file.path(nuccodefileDir,sprintf('nuccode_%s_%gmer_%s.dat',GenomeInfoDb::providerVersion(genome), np, chr));
  nuccodefile = file.path(nuccodefileDir,sprintf('nuccode_%s_%gmer_%s.dat',mm, np, chr));

  if (file.exists(nuccodefile)) {
    #cat(sprintf('readNucleotideCodeForChromosome: loading saved data for %s: refgenome=%s\n', chr, GenomeInfoDb::providerVersion(genome)));
    cat(sprintf('readNucleotideCodeForChromosome: loading saved data for %s: refgenome=%s\n', chr, mm));
    print(nuccodefile);
    nuccode = readRDS(nuccodefile);
  } else {

    if (mapdir != '') {
      unmappable = pickUnmappaleBasesByMappability(mapdir, chr);
    }

    seq=as.character(genome[[chr]])  # convert to a sequence string

    #Splits the sequence into individual characters and initializes an integer vector seqcode of the same length.
    nseq = nchar(seq);
    seqarray=vector("integer", nseq);
    nucleotide= c('A','C','G','T');

    #    cat(sprintf('readNucleotideCodeForChromosome:l79:%s\n', chr));
    sst <- strsplit(seq, "")[[1]]
    #    cat(sprintf('readNucleotideCodeForChromosome:l81:%s\n', chr));
    seqcode=vector(mode="integer",length=nseq)

    #Assigns integer codes to the nucleotide bases in the sequence. Bases that are not 'A', 'C', 'G', or 'T' are set to NA.
    for (ii in 1:length(nucleotide)) {
      seqcode[sst == nucleotide[ii]]=ii;
    }
    seqcode[seqcode==0] = NA;

    #If there are unmappable bases, it masks them in the sequence by setting the corresponding positions in seqcode to NA.
    if (length(unmappable)>1) {
      cat(sprintf('masking  unmappable bases.'));
      lenmap= length(unmappable);
      lenseq= nchar(seq);

      if (lenmap > lenseq) {
        stop('The mappability file is too long.');
      }
      if (lenmap < lenseq) {
        unmappable= c(unmappable, rep(FALSE, length=lenseq-lenmap));
      }
      seqcode[unmappable] = NA;
    }

    #  cat(sprintf('readNucleotideCodeForChromosome:l89:%s\n', chr));
    #np = 4;

    # Calculates the maximum code value and the base values for the nucleotide positions.
    maxcode = 4^np+1;
    bases = 4^((np-1):0);  # 1024,256,64,16,4,1

    #Creates a sequence of indices for shifting the sequence.
    idx=seq(from=-np/2+1, to=np/2);

    #Defines a function shiftsequence that shifts the sequence by a given amount and returns the shifted sequence.
    shiftsequence<-function(shift) {
      rg=(1:nseq)+shift;
      rg[rg<1] = NA;
      rg[rg>nseq] = NA;
      seqcode2 = seqcode[rg];
      seqcode2;
    };

    #Creates a vector codes of the same length as the sequence and initializes it to 0.
    codes = vector(mode="integer",length=nseq);

    #Calculates the code for each position in the sequence.
    for (ii in 1:length(idx)) {
      codes = codes + (shiftsequence(ii)-1) * bases[ii];
    }

    codes = codes + 1;
    codes[is.na(codes)]= maxcode;

    #Calculates the frequency of each nucleotide code and creates a data frame dat with the sequence labels and their counts.
    freq = vector(mode="integer", length=maxcode);
    rl = rle(sort(codes));

    label = getNucleotideString(np);
    label = c(label,'other');
    freq[rl$values] = rl$length;
    dat=data.frame(seq=label, count=freq)

    #Creates a list nuccode containing the nucleotide codes and the frequency table. Saves this list to a file and returns it.
    nuccode=list(code=codes, freqtable = dat);

    saveRDS(nuccode, file = nuccodefile);
  }
  nuccode;
}
