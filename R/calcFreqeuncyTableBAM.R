#' Calculate Frequency Table from BAM File
#'
#' Calculates a frequency table from a BAM file. If a precomputed table exists, it loads it; otherwise, it processes each chromosome, sums counts and reference sequence counts, and saves the result.
#'
#' @param bamfile Character. Path to the BAM file.
#' @param refgenome Character. Reference genome name.
#' @param np Integer. Length of the nucleotide sequence to consider (default: 6).
#' @param mapdir Character. Directory where mappability files are stored.
#' @param shifts Numeric vector. Shifts to apply.
#'
#' @return Data frame. Frequency table with counts and reference sequence counts.
#' @export
#'
#' @examples
#' # calcFreqeuncyTableBAM("sample.bam", "mm10", 6, "mapdir", c(0,0))
calcFreqeuncyTableBAM <- function ( bamfile='',refgenome='', np=6,mapdir='',shifts=shifts) {

  #Creates a unique filename for temporary data storage based on a hash of the input parameters.
  tempdata = sprintf('temp_main_data_%s.dat',digest::digest(list(bamfile, refgenome,np, shifts, mapdir)));

  #Checks if the temporary data file exists. If it does, it loads the frequency table from the file using readRDS. If not, it proceeds to calculate the frequency table
  if (file.exists(tempdata)) {
    freqtable = readRDS(tempdata);
  } else {

    #Initializes genome to NA.
    #Loads chromosome range information from the reference genome using loadChromosomeRange.
    #Extracts the chromosome names from the chromosome information.
    #Loads the reference genome using loadReferenceGenome.
    genome <- NA;
    chrominfo = loadChromosomeRange(refgenome);
    chroms = names(chrominfo);
    genome <- loadReferenceGenome(refgenome);

    # print(chroms);
    #Calculates the frequency table for each chromosome using calcFrequencyTableBAMChromosome and stores the results in freqs.
    #The commented lines show alternative ways to calculate the frequency table, including parallel processing and limiting to specific chromosomes.
    freqs=lapply(chroms, function(x) calcFrequencyTableBAMChromosome(bamfile, x,genome=genome, np=np,mapdir=mapdir,shifts=shifts));
    #freqs=parallel::mclapply(chroms, function(x) calcFrequencyTableBAMChromosome(bamfile, x,genome=genome, np=np,mapdir=mapdir,shifts=shifts), mc.cores=2);
    #freqs=lapply(c('chr1','chr2'), function(x) calcFrequencyTableBAMChromosome(bamfile, x,genome=genome, np=np,mapdir=mapdir,shifts=shifts));

    #Initializes countsum and refseqcountsum with the counts from the first chromosome.
    #If there are multiple chromosomes, it sums the counts and reference sequence counts from all chromosomes.
    countsum = freqs[[1]]$count;
    refseqcountsum = freqs[[1]]$refseqcount;
    if (length(freqs)>1) {
      for (ll in 2:length(freqs)) {
        countsum = countsum + freqs[[ll]]$count;
        refseqcountsum = refseqcountsum + freqs[[ll]]$refseqcount;
      }
    }
    #Creates a data frame freqtable with the summed counts and reference sequence counts.
    #Saves the frequency table to the temporary data file using saveRDS.
    freqtable = data.frame(count=countsum,refseqcount= refseqcountsum, row.names=row.names(freqs[[1]]));
    saveRDS(freqtable, file = tempdata);
  }
  freqtable;
}

#' Load Reference Genome
#'
#' Loads a reference genome as a BSgenome object.
#'
#' @param refgenome Character. Reference genome name ("mm9", "mm10", "hg19","hg38").
#'
#' @return BSgenome object representing the reference genome.
#' @export
#'
#' @examples
#' # loadReferenceGenome("mm10")
loadReferenceGenome <- function(refgenome) {
  if (refgenome == 'mm9') {
    genome <- BSgenome.Mmusculus.UCSC.mm9::BSgenome.Mmusculus.UCSC.mm9;
  } else if (refgenome == 'mm10') {
    genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10;
  } else if (refgenome == 'hg19') {
    genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19;
  } else if (refgenome == 'hg38') {   #
    genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38;
  } else {
    stop(sprintf('Unsupported ref. genome %s\n', refgenome));
  }
  genome;
}



#' Calculate Frequency Table for a Chromosome from BAM File
#'
#' Calculates a frequency table for a single chromosome from a BAM file.
#'
#' @param bamfile Character. Path to the BAM file.
#' @param chr Character. Chromosome name.
#' @param genome BSgenome object. Reference genome.
#' @param np Integer. Length of the nucleotide sequence to consider (default: 6).
#' @param mapdir Character. Directory where mappability files are stored.
#' @param shifts Numeric vector. Shifts to apply (default: c(0,0)).
#'
#' @return Data frame. Frequency table for the chromosome.
#' @export
#'
#' @examples
#' # calcFrequencyTableBAMChromosome("sample.bam", "chr1", genome, 6, "mapdir", c(0,0))

calcFrequencyTableBAMChromosome<-function(bamfile, chr, genome=NA, np=6,mapdir='',shifts=c(0,0)) {

  #tempdata = sprintf('temp_data_%s.dat',digest::digest(list(bamfile, chr, GenomeInfoDb::providerVersion(genome),np, shifts, mapdir)));
  tempdata = sprintf('temp_data_%s.dat',digest::digest(list(bamfile, chr, mm,np, shifts, mapdir)));
  if (file.exists(tempdata)) {
    #  	  cat(sprintf('calcFrequencyTableBAMChromosome: loading saved data for %s: refgenome=%s\n', chr, providerVersion(genome)));
    freqtable = readRDS(tempdata);
  } else {

    #~ 	  if (mappability) {
    #~ 		 mapdir = getMapDirectory(genome);
    #~ 	  } else {
    #~ 		  mapdir ="";
    #~ 	  }

    #unmappable = pickUnmappaleBasesByMappability(mapdir, chr);

    cuts=readCutSitesPerChromFromBAM(bamfile, chr, shifts=shifts);  # ok

    # cat(sprintf('calcFrequencyTableBAMChromosome: finished reading cuts for %s\n', chr));
    nuccode<-readNucleotideCodeForChromosome(chr,nmer=np,genome=genome,mapdir=mapdir); # ok


    cutfreq= rle(cuts);
    tab=cbind(loc=cutfreq$values, count = cutfreq$lengths);
    nucfreq = vector("numeric", length=length(nuccode$freqtable$count));
    codetochange=nuccode$code[cutfreq$values-floor(np/2)];
    for (x in 1:length(codetochange)) {
      nucfreq[codetochange[x]] = nucfreq[codetochange[x]] + cutfreq$lengths[x];
    }

    freqtable = data.frame(count=nucfreq,refseqcount= nuccode$freqtable$count, row.names=nuccode$freqtable$seq);
    cat(sprintf('calcFrequencyTableBAMChromosome: saving the result for %s\n', chr));
    saveRDS(freqtable, file = tempdata);
    # browser()
  }
  freqtable;
}



#' Pick Unmappable Bases by Mappability
#'
#' Identifies unmappable bases in a chromosome using mappability files.
#'
#' @param mapfiledir Character. Directory containing mappability files.
#' @param chr Character. Chromosome identifier.
#' @param map_seqlength Integer. Length of the sequence to map (default: 35).
#'
#' @return Logical vector indicating unmappable bases.
#' @export
#'
#' @examples
#' # pickUnmappaleBasesByMappability("mapdir", "chr1")
pickUnmappaleBasesByMappability<- function(mapfiledir, chr, map_seqlength=35) {

  # Construct the file path for the mappability file
  mapfile = sprintf('%s/%sb.out',  mapfiledir, chr);
  cat(sprintf('reading mappability file: %s...\n', mapfile));

  # Open the mappability file in binary read mode
  fid1 = file(mapfile, 'rb');

  # Get the size of the file
  fileSize <- file.info(mapfile)$size;

  # Read the mappability data from the file
  mappability = readBin(fid1, integer(), n = fileSize, size = 1, signed = FALSE, endian = 'little');

  # Close the file
  close(fid1);

  # Shift the mappability array forward and backward
  forwardcutmappability = shiftarray(mappability, -1);
  backwardcutmappability = shiftarray(mappability, map_seqlength);

  # Determine unmappable bases
  unmappable = (forwardcutmappability != 1) & (backwardcutmappability != 1);

  # Clear variables to free memory
  mappability = {};
  forwardcutmappability = {};
  backwardcutmappability = {};

  # Return the unmappable bases
  unmappable;
}



#' Shift Array
#'
#' Shifts a numeric vector by a specified number of positions. Positions shifted out are replaced with zeros.
#'
#' @param map Numeric vector to be shifted.
#' @param bp Integer. Number of positions to shift (positive: right, negative: left).
#'
#' @return Numeric vector, shifted by bp positions.
#' @export
#'
#' @examples
#' shiftarray(c(1, 2, 3, 4, 5), 2)    # returns c(0, 0, 1, 2, 3)
#' shiftarray(c(1, 2, 3, 4, 5), -2)   # returns c(3, 4, 5, 0, 0)
shiftarray<-function(map, bp) {
  # Convert bp to an integer
  bp = as.integer(bp);
  # Get the length of the map
  ln = length(map);
  # If bp is positive, shift to the right
  if (bp>0) {
    newmap = c(rep(0,bp), map[1:(ln-bp)]);
    # If bp is negative, shift to the left
  } else if (bp<0) {
    newmap = c(map[(-bp+1):ln], rep(0,-bp));
    # If bp is zero, return the original map
  } else {
    newmap = map;
  }
  # Return the new map
  newmap;
}


#' Get Nucleotide String Combinations
#'
#' Generates all possible nucleotide strings of a given length.
#'
#' @param np Integer. Length of the nucleotide strings to generate.
#'
#' @return Character vector of all possible nucleotide strings of length np.
#' @export
#'
#' @examples
#' getNucleotideString(2)
#' getNucleotideString(3)
getNucleotideString<-function(np) {
  # Define the basic nucleotides
  nucleotide = c('A', 'C', 'G', 'T')

  # Initialize the result with the basic nucleotides
  NN = nucleotide

  # If np is greater than or equal to 2, generate combinations
  if (np >= 2) {
    for (i in 2:np) {
      NT = c()
      for (j in 1:length(NN)) {
        # Concatenate each existing string with each nucleotide
        NT = c(NT, paste(NN[j], nucleotide, sep=''))
      }
      # Update NN with the new combinations
      NN = NT
    }
  }

  # Return the final list of nucleotide strings
  NN
}
