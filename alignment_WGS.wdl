#
# alignment_WGS_WF
# Cromwell workflow for aligning WGS data to reference sequence
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

workflow alignment_WGS_WF {
	String genome_ref

	String workDir
	String log
	String tmp

	String inputSamplesFile
	String inputFastqsFile
	Array[Array[String]] inputSamples = read_tsv(inputSamplesFile)
	Array[Array[String]] inputFastqs = read_tsv(inputFastqsFile)

	Int threads
	String bwa
	String bwa_multifastq_script
	String samtools
	String java
	String picard
	String qualimap

	scatter (i in range(length(inputSamples))) {
		call bwaTask {
			input:
				ref = genome_ref,
				name = inputSamples[i][0],
				fastqs = inputFastqs[i],
				outDir = "${workDir}/bam",
				log = log,
				tmp = tmp,
				threads = threads,
				bwa = bwa,
				bwa_multifastq_script = bwa_multifastq_script,
				samtools = samtools,
				java = java,
				picard = picard,
				qualimap = qualimap
		}
	}
}

#######################################################################
task bwaTask {
	String ref
	String name
	Array[String] fastqs
	Int length = length(fastqs)
	String outDir
	String log
	String tmp
	Int threads
	String bwa
	String bwa_multifastq_script
	String samtools
	String java
	String picard
	String qualimap

	Int mem = 2048
	Int mem_qualimap = threads*mem

	String RG = "'@RG\\tID:${name}\\tLB:${name}\\tPL:ILLUMINA\\tPU:${name}\\tSM:${name}'"

	command <<<
		mkdir -p ${outDir} ${log} ${tmp} ;
		BAM="${outDir}/${name}.sorted.bam" ;
		BAMINDEX="${outDir}/${name}.sorted.bai" ;
		if [[ ! -f $BAM || ! -f $BAMINDEX ]]; then
			if [[ ${length} -ge "4" ]]; then
				bash ${bwa_multifastq_script} \
					${name} \
					${outDir} \
					${ref} \
					${bwa} \
					${samtools} \
					${tmp} \
					${threads} \
					${mem} \
					${sep=' ' fastqs} ;
			else
				${bwa} mem \
					-R ${RG} \
					-t ${threads} \
					${ref} ${sep=' ' fastqs} \
					| ${samtools} sort -T ${tmp} -@ ${threads} -m ${mem}M \
					-O BAM -o ${outDir}/${name}.sorted.bam ;
				${samtools} index ${outDir}/${name}.sorted.bam ${outDir}/${name}.sorted.bai ;
			fi ;
		else
			echo -e "${name}.sorted.bam exists" ;
		fi ;

		BAM="${outDir}/${name}.sorted.md.bam" ;
		BAMINDEX="${outDir}/${name}.sorted.md.bai" ;
		if [[ ! -f $BAM || ! -f $BAMINDEX ]]; then
			${java} -Xms${mem_qualimap}M -Xmx${mem_qualimap}M -Djava.io.tmpdir=${tmp} -jar \
				${picard} MarkDuplicates \
				I=${outDir}/${name}.sorted.bam \
				O=${outDir}/${name}.sorted.md.bam \
				M=${outDir}/${name}.metrics.txt \
				CREATE_INDEX=true ;
		else
			echo -e "${name}.sorted.md.bam exists" ;
		fi ;
		
		if ls ${outDir}/${name}_bamqc/${name}_bamqc.pdf 1> /dev/null 2>&1; then
			echo -e "bamqc exists" ;
		else
			${qualimap} bamqc \
				-c \
				-nt ${threads} \
				-outdir ${outDir}/${name}_bamqc \
				-outfile ${name}_bamqc.pdf \
				-outformat PDF:HTML \
				-sd -sdmode 0 \
				-bam ${outDir}/${name}.sorted.md.bam \
				--java-mem-size=${mem_qualimap}M ;
		fi ;

		exit 0 ;
	>>>
	runtime { 
		cpu: "${threads}"
		memory: 2*mem+" MB"
		my_job_name: "BWA_MD_qualimap_${name}"
		my_log: "${log}/BWA_MD_qualimap_${name}"
	}
	output {
		File bam = "${outDir}/${name}.sorted.bam"
		File bamIndex = "${outDir}/${name}.sorted.bai"
		File bamMD = "${outDir}/${name}.sorted.md.bam"
		File bamIndexMD = "${outDir}/${name}.sorted.md.bai"
		File bamqc = "${outDir}/${name}_bamqc/${name}_bamqc.pdf"
	}
}
