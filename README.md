# EVOLINC-I 1.0
-
Evolinc-I is a long intergenic noncoding RNA (lincRNA) identification workflow that also facilitates genome browser visualization of identified lincRNAs and downstream differential gene expression analysis. 

Evolinc-I minimally requires the following input data

1. A set of assembled and merged transcripts from Cuffmerge or Cuffcompare in gene transfer format (GTF)
2. A reference genome (FASTA)
3. A reference genome annotation (GFF)

Optional input data

1. Transposable Elements database (FASTA)
2. Known LincRNA (FASTA)
3. Transcription start site coordinates (BED)
 
# Availablility
-
### Using Docker image

Since there are several dependencies (can be seen in Dockerfile) to Evolinc-I to make it run on your linux or MAC OS, we highly recommend to use Docker image for [Evolinc-I](https://hub.docker.com/r/cyverse/evolinc-i/) or use the [Dockerfile](https://hub.docker.com/r/cyverse/evolinc-i/~/dockerfile/) to build an image and then use the built image for running Evolinc-I

```sudo docker pull cyverse/evolinc-i:1.0```

```sudo docker run --rm -v $(pwd):/working-dir -w /working-dir cyverse/evolinc-i:1.0 -c Sample_cuffcompare_out.gtf -g TAIR10_chr.fasta -r TAIR10_genes.gff -o test_out_new_no_overlap -n 4```

### Using CyVerse Discovery Environment

The [Evolinc-I app](https://de.cyverse.org/de/?type=apps&app-id=e980754e-8050-11e6-97c3-008cfa5ae621&system-id=de) is currently integrated in CyVerse’s Discovery Environment and is free to use by researchers. The complete tutorial is available at this [CyVerse wiki](https://wiki.cyverse.org/wiki/display/TUT/Evolinc+in+the+Discovery+Environment). 

# Issues
-
If you experience any issues with running this Docker image, please contact *upendradevisetty at goooglemail.com* 

# Copyright free
-
The sources in this [Github](https://github.com/Evolinc/Evolinc-I) repository, are copyright free. Thus you are allowed to use these sources in which ever way you like. Here is the full [MIT](https://choosealicense.com/licenses/mit/#) license.

# Citing Evolinc-I
-
Evolinc-I is currently available as a preprint.
