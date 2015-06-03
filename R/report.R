# R program for bam file stats
library(ggplot2)
library(reshape2)
library(grid)

args           <- commandArgs(TRUE)
file_name      <- args[ 1 ]
sample_name    <- args[ 2 ]
png_dir        <- args[ 3 ]
coverage_list  <- args[ 4 ]

if ( FALSE )
{
    file_name   <- '/home/eigos/Tmp/NOG.chiba/summary_R.txt'
    sample_name <- 'NOG samples'
    png_dir     <- '/home/eigos/Tmp'
    coverage_list <- '2,10,20,30,40,50,100'
}

file <- read.table( file_name, header=TRUE, sep="\t" )

#
# Draw
#

#
# Mean insert
#
gg <- ggplot( file, aes( bam_filename, median_insert_size, fill=bam_filename) )
gg <- gg + geom_bar( width=0.7, stat="identity" )
gg <- gg + geom_point( size = 5, shape = 8, aes( bam_filename, mean_insert_size) )
gg <- gg + geom_errorbar( width=0.3, aes(ymin=median_insert_size - insert_size_sd, ymax=median_insert_size + insert_size_sd ) )
gg <- gg + ggtitle( sample_name )
gg <- gg + theme(plot.title = element_text(size = 30),
                 axis.title=element_text(size = 30),
                 axis.text.x=element_text(size = 30, angle = 70, hjust = 1),
                 axis.text.y=element_text(size = 30),
                 legend.text=element_text(size = 30, hjust = 1),
                 legend.key.size = unit( 1, "cm" ),
                 legend.title=element_text(size = 30, hjust = 0))
gg <- gg + scale_fill_hue( l = 45 )
ggsave(file=paste(png_dir, '/', sample_name, '_mean_insert.png', sep=""), plot=gg, dpi=100, width=24,height=18)

#
# Mean depth
#
gg <- ggplot( file, aes( x = bam_filename, y = average_depth, fill=bam_filename) )
gg <- gg + geom_bar( width=0.7, stat="identity" )
gg <- gg + geom_errorbar( data = file, aes( x = bam_filename, ymin = average_depth, ymax = average_depth + depth_stdev ), width = 0.4 )
gg <- gg + ggtitle( sample_name )
gg <- gg + theme(plot.title = element_text(size = 30),
                 axis.title=element_text(size = 30),
                 axis.text.x=element_text(size = 30, angle = 70, hjust = 1),
                 axis.text.y=element_text(size = 30),
                 legend.text=element_text(size = 30, hjust = 1),
                 legend.key.size = unit( 1, "cm" ),
                 legend.title=element_text(size = 30, hjust = 0))
gg <- gg + scale_fill_hue( l = 45 )
ggsave(file=paste(png_dir, '/', sample_name, '_depth.png', sep=""), plot=gg, dpi=100, width=24,height=18)

#
# coverage
#
cov_list = strsplit( coverage_list, ',' )
data_list = c( 'bam_filename' )
for ( i in (1:length( cov_list[[1]] ) ) ) {
    cov_str <- paste('X', cov_list[[1]][i], 'x_ratio', sep='' )
    data_list <- c( data_list, c = cov_str)
}
df = melt( file[ data_list ] )

gg <- ggplot( )
gg <- gg + ylim( 0, 1 )
gg <- gg + geom_line( data = df, width=0.7, aes( bam_filename, value, group = variable, colour = variable ) )
gg <- gg + geom_point( data = df, size=5, aes( bam_filename, value, group = variable, colour = variable ) )
gg <- gg + ylab( "Coverage" )
gg <- gg + ggtitle( sample_name )
gg <- gg + theme(plot.title = element_text(size = 30),
                 axis.title=element_text(size = 30),
                 axis.text.x=element_text(size = 30, angle = 70, hjust = 1),
                 axis.text.y=element_text(size = 30),
                 legend.text=element_text(size = 30, hjust = 1),
                 legend.key.size = unit( 1, "cm" ),
                 legend.title=element_text( size = 30, hjust = 0) )
gg <- gg + scale_color_discrete(name="Coverage" )
gg <- gg + scale_fill_hue( l = 45 )
ggsave(file=paste(png_dir, '/', sample_name, '_coverage_line.png', sep=""), plot=gg, dpi=100, width=24,height=18)

gg <- ggplot( )
gg <- gg + ylim( 0, 1 )
gg <- gg + geom_bar( data = df, width=0.7, stat="identity", position="dodge", aes( bam_filename, value, group = variable, fill = variable ) )
gg <- gg + ylab( "Coverage" )
gg <- gg + ggtitle( sample_name )
gg <- gg + theme(plot.title = element_text(size = 30),
                 axis.title=element_text(size = 30),
                 axis.text.x=element_text(size = 30, angle = 70, hjust = 1),
                 axis.text.y=element_text(size = 30),
                 legend.text=element_text(size = 30, hjust = 1),
                 legend.key.size = unit( 1, "cm" ),
                 legend.title=element_text( size = 30, hjust = 0) )
gg <- gg + scale_color_discrete(name="Coverage" ) 
gg <- gg + scale_fill_hue( l = 45 )
ggsave(file=paste(png_dir, '/', sample_name, '_coverage_dodge.png', sep=""), plot=gg, dpi=100, width=24,height=18)

data_list = c( 'bam_filename' )
for ( i in (1:length( cov_list[[1]] ) ) ) {
    new_str <- paste( 'x', cov_list[[1]][ i ], sep='' )
    old_A_str <- paste( 'X', cov_list[[1]][ i ], 'x_ratio', sep='' )
    print( i )
    print( new_str )
    print( old_A_str )
    if ( i < length( cov_list[[1]] ) ) {
        old_B_str <- paste( 'X', cov_list[[1]][ i + 1 ], 'x_ratio', sep='' )
        print( old_B_str )
        file[ new_str ] = file[ old_A_str ] - file[ old_B_str ]
    } else {
        file[ new_str ] = file[ old_A_str ]
    }

    data_list <- c( data_list, c = new_str )
}
cov <- melt( file[ data_list ] )
gg <- ggplot( cov, aes( bam_filename, value, group = variable, fill = variable )  )
gg <- gg + geom_bar( position="fill", width=0.7, stat="identity" )
gg <- gg + ggtitle( sample_name )
gg <- gg + theme(plot.title = element_text(size = 30),
                 axis.title=element_text(size = 30),
                 axis.text.x=element_text(size = 30, angle = 70, hjust = 1),
                 axis.text.y=element_text(size = 30),
                 legend.text=element_text(size = 30, hjust = 1),
                 legend.key.size = unit( 1, "cm" ),
                 legend.title=element_text(size = 30, hjust = 0))
gg <- gg + scale_fill_hue( l = 45 )
ggsave(file=paste(png_dir, '/', sample_name, '_coverage_bar.png', sep=""), plot=gg, dpi=100, width=24,height=18)

#
# Duplicates
#
file$duplicate_ratio = file$num_duplicate_reads / file$num_total_reads
gg <- ggplot( file, aes( x = bam_filename, y = duplicate_ratio, fill=bam_filename) )
gg <- gg + geom_bar( width=0.7, stat="identity" )
gg <- gg + ggtitle( sample_name )
gg <- gg + theme(plot.title = element_text(size = 30),
                 axis.title=element_text(size = 30),
                 axis.text.x=element_text(size = 30, angle = 70, hjust = 1),
                 axis.text.y=element_text(size = 30),
                 legend.text=element_text(size = 30, hjust = 1),
                 legend.key.size = unit( 1, "cm" ),
                 legend.title=element_text(size = 30, hjust = 0))
gg <- gg + scale_fill_hue( l = 45 )
ggsave(file=paste(png_dir, '/', sample_name, '_duplicate.png', sep=""), plot=gg, dpi=100, width=24,height=18)


#
# Alignment
#
gg <- ggplot( file, aes( bam_filename, num_total_reads, fill = bam_filename ) )
gg <- gg + geom_bar( width=0.7, stat="identity" )
gg <- gg + ggtitle( sample_name )
gg <- gg + theme(plot.title = element_text(size = 30),
                 axis.title=element_text(size = 30),
                 axis.text.x=element_text(size = 30, angle = 70, hjust = 1),
                 axis.text.y=element_text(size = 30),
                 legend.text=element_text(size = 30, hjust = 1),
                 legend.key.size = unit( 1, "cm" ),
                 legend.title=element_text(size = 30, hjust = 0))
gg <- gg + scale_fill_hue( l = 45 )
ggsave(file=paste(png_dir, '/', sample_name, '_metrics.png', sep=""), plot=gg, dpi=100, width=24,height=18)
