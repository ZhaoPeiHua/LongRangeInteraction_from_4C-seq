# LongRangeInteraction_from_4C-seq  
The R script (LongRangeInteraction.R) for detecting long-range interactions (far-cis and trans interactions) from 4C-seq.  
  
# **PREREQUISITE** :  
R software (any version).  
      
# **ARGUMENTS** :  
1. filename : A tab-delimited text file containing the number of mapped 4C-seq reads per fragment end. The file has to be sorted by chromosome and then by position of the first restriction enzyme sites in ascending order. Here is an example of the text file.  
   * 1st column : Chromosome.    
   * 2nd column : Position of the first RE sites in the genome.    
   * 3rd column : Number of reads that start at the first RE sites.

    |  Chromosome   |  Position   | Reads |
    |-----|-----|-----|
    |chr1  |11159   |0|
    |chr1  |12410   |0|
    |chr1  |12460   |0|
    |...  |...   |...|
    |chr1  |1457897   |0|
    |chr1  |1458262   |12|
    |chr1  |1458337   |0|

    
    
1. vp.chr : The chromosome where the viewpoint locates.  
1. vp.pos : The position of the viewpoint, could be any coordinate inside the viewpoint.  
      
# **OPTIONS** :  
1. window_size : The number of fragment ends to be defined as an interaction. Default=400.  
1. min_reads : The minimum number of reads is required for fragment end to be defined as an observed fragment. Default=2.  
1. masked_region_length : Near-cis region (+/- 'masked_region_length' from the viewpoint) will be masked out. Default=5000000.  
      
      
# **USAGE** :  
Open R software  
source("/path/to/LongRangeInteraction.R")  
example_data="/path/to/example_data.txt"  
LongRangeInteraction(filename=example_data,vp.chr="chr9",vp.pos=21969726,window_size=400,min_reads=2,masked_region_length=5000000)  
