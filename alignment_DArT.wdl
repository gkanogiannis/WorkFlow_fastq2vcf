#
# alignment_DArT_WF
# Cromwell workflow for aligning DArT data to reference sequence
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

workflow alignment_DArT_WF {
	String genome_ref

	String workDir
	String log
	String tmp

	String inputSamplesFile
	Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)
	
	String threads
	String bwa
	String samtools

	scatter (sample in inputSamples) {
		String name = sample[0]
		
		call bwaMem {
			input:
				ref = genome_ref,
				name = name,
				fastq1 = sample[1],
				outDir = workDir,
				log = log,
				tmp = tmp,
				threads = threads,
				bwa = bwa,
				samtools = samtools
		}
	}
}

task bwaMem {
	String ref
	String name
	File fastq1
	String outDir
	String log
	String tmp
	String threads
	String bwa
	String samtools

	String RG = "'@RG\\tID:${name}\\tLB:${name}\\tPL:ILLUMINA\\tPU:${name}\\tSM:${name}'"

	command <<<
		mkdir -p ${outDir} ${tmp} ;
		BAM="${outDir}/${name}.sorted.bam" ;
		BAMINDEX="${outDir}/${name}.sorted.bai" ;
		if [[ ! -f $BAM || ! -f $BAMINDEX ]]; then
			${bwa} mem \
				-R ${RG} \
				-t ${threads} \
				${ref} ${fastq1} \
				| ${samtools} sort -T ${tmp} -@ ${threads} \
				-O BAM -o ${outDir}/${name}.sorted.bam ;
			${samtools} index ${outDir}/${name}.sorted.bam ${outDir}/${name}.sorted.bai ;
		fi ;
	>>>
	runtime { 
		cpu: "${threads}"
		memory: "16000 MB"
		my_job_name: "BWA_${name}"
		my_log: "${log}/BWA_${name}"
	}
	output {
		File bam = "${outDir}/${name}.sorted.bam"
		File bamIndex = "${outDir}/${name}.sorted.bam.bai"
	}
}
