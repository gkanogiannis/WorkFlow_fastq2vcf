[![DOI](https://zenodo.org/badge/433710013.svg)](https://zenodo.org/badge/latestdoi/433710013)

# WorkFlow_fastq2vcf
 A Cromwell workflow to call and filter variants from sequencing data.

1. QC on samples raw sequencing data
````
java -Dconfig.file=sge.config -Xms32G -Xmx32G -jar cromwell.jar run fastQC.wdl --inputs fastQC.inputs.json > fastQC.log
````
2. Align samples raw sequencing data (fastq -> bam)
````
java -Dconfig.file=sge.config -Xms32G -Xmx32G -jar cromwell.jar run alignment_WGS.wdl --inputs alignment_WGS.inputs.json > alignment_WGS.log
````
3. Call variants (bam -> vcf)
```
java -Dconfig.file=sge.config -Xms32G -Xmx32G -jar cromwell.jar run VariantDiscovery.wdl --inputs VariantDiscovery.inputs.json > VariantDiscovery.log
```
4. Filter variants
````
java -Dconfig.file=sge.config -Xms8G -Xmx8G -jar cromwell.jar run VariantFiltering.wdl --inputs VariantFiltering.inputs.json > VariantFiltering.log
````
