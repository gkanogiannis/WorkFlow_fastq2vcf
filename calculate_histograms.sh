#
# calculate_stats.sh
# Bash script for calculating histograms of Ho, He, PIC, MAF, MissPercentage from (g)vcf
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

#$1 input stats
#$2 output basename

unset R_HOME
Rscript --vanilla <(echo -e "\
library(reshape2)
library(ggplot2)
library(lattice)

data <- read.table(file=gzfile('$1'),header=T,stringsAsFactors=F,fill=T)
gg <- melt(data, id.vars=c('CHR','POS'), measure.vars=c('Ho','He','PIC','MAF','MissPercentage'), variable.name='Values')

pdf(file=paste0('$2','.Ho_He_PIC_MAF_Mis.hist.pdf'),paper='usr',width=11, height=8.5)
print(ggplot(gg, aes(x=value, fill=Values)) +
    geom_histogram(binwidth=0.005)+
    facet_grid(Values~.) + 
    scale_y_continuous(expand=c(0, 0)) +
    scale_x_continuous(expand=c(0, 0), limits=c(0,1)))
  dev.off()
")
