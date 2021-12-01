#
# HaplotypeCaller_WF
# Cromwell workflow for running HaplotypeCaller of GATK
#
# Copyright (C) 2021 Anestis Gkanogiannis <anestis@gkanogiannis.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
#

workflow HaplotypeCaller_WF {
	String genome_ref

	String name
	File bam
	File bamIndex

	String listChromosomeFile
	Array[String] listChromosome = read_lines(listChromosomeFile)
	String? listScaffoldFile
	Boolean withScaffolds

	String basename

	String workDir
	String tmp
	String log

	String threads
	String gatk

	scatter (chromosome in listChromosome) {
		call HaplotypeCaller as HCChromosome {
			input: 
				ref = genome_ref,
				name = name,
				bam = bam,
				bamIndex = bamIndex,
				outDir = workDir,
				basename=basename,
				log = log,
				tmp = tmp,
				interval = chromosome,
				suffix = chromosome,
				threads = 2,
				gatk = gatk
		}
	}

	if(withScaffolds){
		call HaplotypeCaller as HCScaffold {
			input: 
				ref = genome_ref,
				name = name,
				bam = bam,
				bamIndex = bamIndex,
				outDir = workDir,
				basename=basename,
				log = log,
				tmp = tmp,
				interval = listScaffoldFile,
				suffix = "scaffold",
				threads = 2,
				gatk = gatk
		}
	}

	if(withScaffolds){
		call CombineGVCFs_withScaffolds {
			input:
				ref = genome_ref,
				name = name,
				outDir = workDir,
				basename=basename,
				log = log,
				tmp = tmp,
				gvcfsChromosome = HCChromosome.gvcf,
				gvcfsIndexChromosome = HCChromosome.gvcfIndex,
				gvcfScaffold = HCScaffold.gvcf,
				gvcfIndexScaffold = HCScaffold.gvcfIndex,
				threads = 2,
				gatk = gatk
		}
	}
	if(!withScaffolds){
		call CombineGVCFs_withoutScaffolds {
			input:
				ref = genome_ref,
				name = name,
				outDir = workDir,
				basename=basename,
				log = log,
				tmp = tmp,
				gvcfsChromosome = HCChromosome.gvcf,
				gvcfsIndexChromosome = HCChromosome.gvcfIndex,
				threads = 2,
				gatk = gatk
		}
	}

	output {
		File gvcf = "${workDir}/${name}.g.vcf.gz"
		File gvcfIndex = "${workDir}/${name}.g.vcf.gz.tbi"
	}
}

#################################################################
task HaplotypeCaller {
	String ref
	String name
	File bam
	File bamIndex
	String outDir
	String basename
	String log
	String tmp
	String interval
	String suffix
	String gatk

	Int? threads=2
	Int mem=7000
	Int mems=3500
	
	# --min-base-quality-score 28 \
	command <<<
		mkdir -p ${outDir} ${tmp} ;
		GVCF="${outDir}/${name}_${suffix}.g.vcf.gz" ;
		GVCFINDEX="${outDir}/${name}_${suffix}.g.vcf.gz.tbi" ;
		if [[ ! -f $GVCF || ! -f $GVCFINDEX ]]; then
			${gatk} --java-options "-Xms${mems}M -Xmx${mem}M" HaplotypeCaller \
				-R ${ref} \
				-I ${bam} \
				-O ${outDir}/${name}_${suffix}.g.vcf.gz \
				-ERC GVCF \
				-G StandardAnnotation \
				-G StandardHCAnnotation \
				--do-not-run-physical-phasing \
				-L ${interval} \
				--min-base-quality-score 28 \
				--kmer-size 11 --kmer-size 17 --kmer-size 23 \
				--kmer-size 29 --kmer-size 35 --kmer-size 41 ;
		else
			echo -e "${name}_${suffix}.g.vcf.gz exists" ;
		fi ;
		exit 0 ;
	>>>
	runtime { 
		cpu: threads
		memory: 2.0*mem/threads + " MB"
		my_job_name: "HC_${basename}_${name}_${suffix}"
		my_log: "${log}/HC_${basename}_${name}_${suffix}"
		maxRetries: 1
	}
	output {
		File gvcf = "${outDir}/${name}_${suffix}.g.vcf.gz"
		File gvcfIndex = "${outDir}/${name}_${suffix}.g.vcf.gz.tbi"
	}
}

task CombineGVCFs_withScaffolds {
	String ref
	String name
	String outDir
	String basename
	String log
	String tmp
	Array[File] gvcfsChromosome
	Array[File] gvcfsIndexChromosome
	File gvcfScaffold
	File gvcfIndexScaffold
	String gatk

	Int threads
	Int mem=4000
	Int mems=2500

	command <<<
		mkdir -p ${outDir} ${tmp} ;
		GVCF="${outDir}/${name}.g.vcf.gz" ;
		GVCFINDEX="${outDir}/${name}.g.vcf.gz.tbi" ;
		if [[ ! -f $GVCF || ! -f $GVCFINDEX ]]; then
			${gatk} --java-options "-Xms${mems}M -Xmx${mem}M" CombineGVCFs \
				-R ${ref} \
				-O ${outDir}/${name}.g.vcf.gz \
				-G StandardAnnotation \
				-G StandardHCAnnotation \
				--variant ${sep=' --variant ' gvcfsChromosome} \
				--variant ${gvcfScaffold} ;
		else
			echo -e "${name}.g.vcf.gz exists" ;
		fi ;
		if [[ -f $GVCF && -f $GVCFINDEX ]]; then
			rm -f `readlink ${gvcfScaffold}` `readlink ${gvcfIndexScaffold}` ;
			rm -f `readlink ${sep='` ; rm -f `readlink ' gvcfsChromosome}` ;
			rm -f `readlink ${sep='` ; rm -f `readlink ' gvcfsIndexChromosome}` ;
			rm -f ${gvcfScaffold} ${gvcfIndexScaffold} ;
			rm -f ${sep=' ; rm -f ' gvcfsChromosome} ;
			rm -f ${sep=' ; rm -f ' gvcfsIndexChromosome} ;
		fi ;
		exit 0 ;
	>>>
	runtime { 
		cpu: threads
		memory: 2.0*mem/threads + " MB"
		my_job_name: "HC_${basename}_CombineGVCFs_${name}"
		my_log: "${log}/HC_${basename}_CombineGVCFs_${name}"
	}
	 output {
		File gvcf = "${outDir}/${name}.g.vcf.gz"
		File gvcfIndex = "${outDir}/${name}.g.vcf.gz.tbi"
	}
}

task CombineGVCFs_withoutScaffolds {
	String ref
	String name
	String outDir
	String basename
	String log
	String tmp
	Array[File] gvcfsChromosome
	Array[File] gvcfsIndexChromosome
	String gatk

	Int threads
	Int mem=4000
	Int mems=2500

	command <<<
		mkdir -p ${outDir} ${tmp} ;
		GVCF="${outDir}/${name}.g.vcf.gz" ;
		GVCFINDEX="${outDir}/${name}.g.vcf.gz.tbi" ;
		if [[ ! -f $GVCF || ! -f $GVCFINDEX ]]; then
			${gatk} --java-options "-Xms${mems}M -Xmx${mem}M" CombineGVCFs \
				-R ${ref} \
				-O ${outDir}/${name}.g.vcf.gz \
				-G StandardAnnotation \
				-G StandardHCAnnotation \
				--variant ${sep=' --variant ' gvcfsChromosome} ;
		else
			echo -e "${name}.g.vcf.gz exists" ;
		fi ;
		if [[ -f $GVCF && -f $GVCFINDEX ]]; then
			rm -f `readlink ${sep='` ; rm -f `readlink ' gvcfsChromosome}` ;
			rm -f `readlink ${sep='` ; rm -f `readlink ' gvcfsIndexChromosome}` ;
			rm -f ${sep=' ; rm -f ' gvcfsChromosome} ;
			rm -f ${sep=' ; rm -f ' gvcfsIndexChromosome} ;
		fi ;
		exit 0 ;
	>>>
	runtime { 
		cpu: threads
		memory: 2.0*mem/threads + " MB"
		my_job_name: "HC_${basename}_CombineGVCFs_${name}"
		my_log: "${log}/HC_${basename}_CombineGVCFs_${name}"
	}
	 output {
		File gvcf = "${outDir}/${name}.g.vcf.gz"
		File gvcfIndex = "${outDir}/${name}.g.vcf.gz.tbi"
	}
}
