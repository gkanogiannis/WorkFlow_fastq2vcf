#!/bin/bash

#
# bwa_multifastq_script.sh
# Bash script for aligning with BWA, multi-fastq samples to reference sequence
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

arr=( "$@" ) ;

name="${arr[0]}" ;
outDir="${arr[1]}" ;
ref="${arr[2]}" ;
bwa="${arr[3]}" ;
samtools="${arr[4]}" ;
tmp="${arr[5]}" ;
threads="${arr[6]}" ;
mem="${arr[7]}" ;

RG="@RG\\tID:${name}\\tLB:${name}\\tPL:ILLUMINA\\tPU:${name}\\tSM:${name}" ;
echo -e "$RG";

length=${#arr[@]} ;
length=$((length-8)) ;

times=$(($length / 2)) ;
echo -e "greater\ttimes=$times" ;
for k in `seq 1 $times`; do
	i=$((2*$k-2));
	j=$((2*$k-1));
	echo -e "$k\t$i\t$j";
	i=$((i+8)) ;
	j=$((j+8)) ;
	${bwa} mem \
		-R ${RG} \
		-t ${threads} \
		${ref} ${arr[i]} ${arr[j]} \
		| ${samtools} sort -T ${tmp} -@ ${threads} -m ${mem}M \
		-O BAM -o ${outDir}/${name}.part_${k}.sorted.bam ;
done ;

${samtools} merge ${outDir}/${name}.sorted.bam ${outDir}/${name}.part_*.sorted.bam ;
${samtools} index ${outDir}/${name}.sorted.bam ${outDir}/${name}.sorted.bai ;
rm -rf ${outDir}/${name}.part_*.sorted.bam ;
