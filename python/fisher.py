#!/usr/local/package/python2.7/2.7.2/bin/python
"""

bamfilter.py


"""
import sys
import os
import re
import pysam
import scipy.special
from scipy.stats import fisher_exact as fisher
import argparse
import logging

#
# Globals
#
f_print = True

arg = None
target = None
remove_chr = None
filter_quals = None

#
# Class definitions
#

############################################################
class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
        return value


#
# Subroutines
#

# data_pair IDs
POS_CHR = 0
POS_COORD = 1
POS_REF = 2
POS_DATA1 = 3
POS_DATA2 = 4
POS_DATA1_FISHER_INS = 5
POS_DATA1_FISHER_DEL = 6
POS_DATA2_FISHER_INS = 7
POS_DATA2_FISHER_DEL = 8
POS_FISHER_SNV = 9
POS_COUNT = 10

############################################################
def Pileup_out( mpileup, w, threshold, mismatch_rate, min_depth ):

    #
    # mpileup format
    #
    # chr1 272 T 24  ,.$.....,,.,.,...,,,.,..^+. <<<+;<<<<<<<<<<<=<;<;7<&
    #
    # 0 chromosome,
    # 1 1-based coordinate,
    # 2 reference base,
    # 3 the number of reads covering the site (1)
    # 4 read bases (1)
    # 5 base qualities (1)
    # 6 the number of reads covering the site (2)
    # 7 read bases (2)
    # 8 base qualities (2)
    #
    global target
    global remove_chr
    global filter_quals

    try:
        #
        # Prepare mpileup data
        #
        mp_list = str( mpileup.translate( None, '\n' ) ).split( '\t' )
        mp_list_len = len( mp_list )
        ref_base_U = mp_list[ 2 ].upper()
        coordinate = mp_list[ 0:3 ]

        #
        # skip if depth is 0
        #
        if mp_list[ 3 ] == '0' or ( mp_list_len > 6 and mp_list[ 6 ] == '0' ):
            return None

        #
        # skip if reference is 'N'
        #
        if ref_base_U == 'N':
            return None

        #
        # data_pair IDs
        # POS_CHR = 0
        # POS_COORD = 1
        # POS_REF = 2
        # POS_DATA1 = 3
        # POS_DATA2 = 4
        # POS_FISHER = 5
        # POS_COUNT = 6
        #
        data_pair = [ mp_list[ 0 ],
                      mp_list[ 1 ],
                      mp_list[ 2 ],
                      AutoVivification(),
                      AutoVivification(),
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0]


        #
        # Loop for 2 bam file case
        #
        if mp_list_len > 6:
            comparison = True
        else:
            comparison = False

        if comparison:
            input_list = [ ( POS_DATA1, mp_list[ 3 ], mp_list[ 4 ], mp_list[ 5 ] ),
                           ( POS_DATA2, mp_list[ 6 ], mp_list[ 7 ], mp_list[ 8 ] ) ]
        else:
            input_list = [ ( POS_DATA1, mp_list[ 3 ], mp_list[ 4 ], mp_list[ 5 ] ) ]

        for data_id, depth, read_bases, qual_list in input_list:

            indel = AutoVivification()

            #
            # Look for deletion/insertion and save info in 'indel' dictionary
            #
            #   ([\+\-])[0-9]+[ACGTNacgtn]+
            #
            # m.group(1): + or - (deletion/insertion)
            # m.group(2): number of deletion/insertion
            # m.group(3): nucleotides
            #
            deleted = 0
            iter = target.finditer( read_bases )
            for m in iter:
                site = m.start()
                type = m.group( 1 )
                num = m.group( 2 )
                bases = m.group( 3 )[ 0:int( num ) ]
                if bases.islower():
                    strand = ( '-1', '1' )
                else:
                    strand = ( '1', '-1' )

                key = '\t'.join( coordinate + [ bases.upper() ] )
                if key in indel[ type ]:
                    indel[ type ][ key ][ strand[ 0 ] ] += 1
                else:
                    indel[ type ][ key ][ strand[ 0 ] ] = 1
                    indel[ type ][ key ][ strand[ 1 ] ] = 0

                read_bases = read_bases[ 0:site - deleted ] + read_bases[ site + int(num) + len( num ) + 1 - deleted: ]
                deleted += 1 + len( num ) + int( num )

            #
            # Remove '^.' and '$'
            #
            read_bases = remove_chr.sub( '', read_bases )
            read_bases = read_bases.translate( None, '$' ) 

            #
            # Error check
            #
            if len( read_bases ) != len( qual_list ):
                logging.error( "mpileup data is not good: {0}, {1}".format( mpileup, read_bases ) )
                return None

            #
            # Count mismatch
            #
            if int( depth ) > min_depth:

                read_bases = read_bases.replace( '.', ref_base_U )
                read_bases = read_bases.replace( ',', ref_base_U.lower() )

                base_num = {
                    "total_A": 0,
                    "total_C": 0,
                    "total_G": 0,
                    "total_T": 0,
                    "A": 0,
                    "C": 0,
                    "G": 0,
                    "T": 0,
                    "a": 0,
                    "c": 0,
                    "g": 0,
                    "t": 0
                }

                #
                # Count number
                #
                for nuc, qual in zip( read_bases, qual_list ):
                    if qual in filter_quals:
                        if nuc in 'ATGCatgc':
                            base_num[ nuc ] += 1
                            base_num[ 'total_' + nuc.upper() ] += 1

                data_pair[ data_id ][ 'bases' ] = read_bases
                data_pair[ data_id ][ 'depth' ] = depth
                ref_num = base_num[ 'total_' + ref_base_U ]

                mis_num = 0
                mis_base_U = None
                for nuc in ( 'A', 'C', 'G', 'T' ):
                    data_pair[ data_id ][ nuc ] = base_num[ nuc ]
                    tmp = nuc.lower()
                    data_pair[ data_id ][ tmp ] = base_num[ tmp ]
                    tmp = 'total_' + nuc
                    data_pair[ data_id ][ tmp ] = base_num[ tmp ]

                    if nuc != ref_base_U:
                        if base_num[ tmp ] > mis_num:
                            mis_num = base_num[ tmp ]
                            mis_base_U = nuc

                #
                # count ins and del
                #
                for type in ( '+', '-' ):
                    if type in indel:
                        for key in indel[ type ].keys():
                            bases = key.split( '\t' )[ 3 ]
                            data_pair[ data_id ][ 'indel' ][ type ][ bases ][ '1' ] = indel[ type ][ key ][ '1' ]
                            data_pair[ data_id ][ 'indel' ][ type ][ bases ][ '-1' ] = indel[ type ][ key ][ '-1' ]
                            data_pair[ data_id ][ 'indel' ][ type ][ bases ][ '0' ] = ( indel[ type ][ key ][ '-1' ] +
                                                                                        indel[ type ][ key ][ '1' ] )


                #
                # Calculate ratio
                #
                data_pair[ data_id ][ 'mis_rate' ] = mis_num / float( depth )
                data_pair[ data_id ][ 'mis_base' ] = mis_base_U
                if mis_base_U:
                    data_pair[ data_id ][ 's_ratio' ]  = float( base_num[ mis_base_U ] ) / ( base_num[ mis_base_U ] + base_num[ mis_base_U.lower() ] )
                else:
                    data_pair[ data_id ][ 's_ratio' ]  = 0

                #
                # Beta distribution for SNV
                #
                data_pair[ data_id ][ '0.1' ] = scipy.special.btdtri( mis_num + 1, ref_num + 1, 0.1 )
                data_pair[ data_id ][ 'mid' ] = ( mis_num + 1 ) / float(ref_num + mis_num + 2 )
                data_pair[ data_id ][ '0.9' ] = scipy.special.btdtri( mis_num + 1, ref_num + 1, 0.9 )
                data_pair[ POS_COUNT ] += 1


        #
        # Fisher
        #
        data_pair[ POS_FISHER_SNV ] = '-'
        data_pair[ POS_DATA1_FISHER_INS ] = '-'

        if data_pair[ POS_COUNT ] == 2 and ref_base_U and mis_base_U:
            # SNV
            odds_ratio, fisher_pvalue = fisher(
                        ( ( int( data_pair[ POS_DATA2 ][ 'total_' + ref_base_U ] ),
                            int( data_pair[ POS_DATA1 ][ 'total_' + ref_base_U ] ) ),
                          ( int( data_pair[ POS_DATA2 ][ 'total_' + mis_base_U ] ),
                            int( data_pair[ POS_DATA1 ][ 'total_' + mis_base_U ] ) ) ),
                        alternative='two-sided'
                )
            data_pair[ POS_FISHER_SNV ] = fisher_pvalue


    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        logging.error( ("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) ) )

    return data_pair


############################################################
def print_data( data, w, min_depth, mismatch_rate ):

    #print data[1]
    # chr pos ref
    str_list = [ '\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-',
                 '\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-' ]
    str_id = 0
    f_print = False
    f_print_indel = False

    for data_id in range( POS_DATA1, POS_DATA1 + data[ POS_COUNT ] ):
        if not data[ data_id ].has_key( 'bases' ):
            continue

        indel = [ '', '', '', '' ]
        for type in data[ data_id ][ 'indel' ].keys():
            if type == '+':
                indel_id = 0
            elif type == '-':
                indel_id = 2

            for bases in data[ data_id ][ 'indel' ][ type ].keys():
                for strand in data[ data_id ][ 'indel' ][ type ][ bases ].keys():
                    if strand == '-1':
                        indel_id_tmp = indel_id + 1
                    elif strand == '1':
                        indel_id_tmp = indel_id
                    elif strand == '0':
                        continue

                    if data[ data_id ][ 'indel' ][ type ][ bases ][ strand ] > 0:
                        f_print_indel = True
                        indel[ indel_id_tmp ] += "{0}:{1},".format(
                                bases,
                                data[ data_id ][ 'indel' ][ type ][ bases ][ strand ] )


        # obs depth A a C c G g T t ins  del  A C G T mis s_ratio 0.1 ratio 0.9
        str_list[ str_id ] = '\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\t{21}\t{22}'.format(
                    data[ data_id ][ 'mis_base' ],
                    data[ data_id ][ 'depth' ],
                    data[ data_id ][ 'A' ],
                    data[ data_id ][ 'a' ],
                    data[ data_id ][ 'C' ],
                    data[ data_id ][ 'c' ],
                    data[ data_id ][ 'G' ],
                    data[ data_id ][ 'g' ],
                    data[ data_id ][ 'T' ],
                    data[ data_id ][ 't' ],
                    indel[ 0 ],
                    indel[ 1 ],
                    indel[ 2 ],
                    indel[ 3 ],
                    data[ data_id ][ 'total_A' ],
                    data[ data_id ][ 'total_C' ],
                    data[ data_id ][ 'total_G' ],
                    data[ data_id ][ 'total_T' ],
                    data[ data_id ][ 'mis_rate' ],
                    data[ data_id ][ 's_ratio' ],
                    data[ data_id ][ '0.1' ],
                    data[ data_id ][ 'mid' ],
                    data[ data_id ][ '0.9' ]
                    )
        str_id += 1

        if (  ( data[ data_id ][ 'depth' ]      >   min_depth       and
                data[ data_id ][ 'mis_base' ]   !=  'N'             and
                data[ data_id ][ 'mis_base' ]   !=  None            and
                data[ data_id ][ 'mis_base' ]   !=  data[ POS_REF ] and
                data[ data_id ][ 'mis_rate']    >   mismatch_rate       )
                or
                f_print_indel
           ):
            f_print = True or f_print

    if f_print:
        outstr = '\t'.join( data[POS_CHR:POS_DATA1] ) + str_list[ 0 ]
        if data[ POS_COUNT ] == 2:
            outstr +=  str_list[ 1 ] + '\t' + str( data[ POS_FISHER_SNV ] )
        outstr +=  '\n'

        w.write( outstr )


############################################################
def Pileup_and_count(
        in_bam1 = None,
        in_bam2 = None,
        out_file = None,
        input_mpileup = None,
        ref_fa = None,
        threshold = 15,
        mismatch_rate = 0.07,
        min_depth = 9
        ):

    global arg
    global f_print
    global target
    global remove_chr
    global filter_quals

    try:

        #
        # Initalize filter quality values
        #
        filter_quals = ''
        for qual in range( 33 + threshold, 33 + threshold + 50 ):
            filter_quals += str( unichr( qual ) )

        #
        # Setup regular expression
        # ([\+\-])[0-9]+[ACGTNacgtn]+
        #
        target = re.compile( '([\+\-])([0-9]+)([ACGTNacgtn]+)' )
        remove_chr = re.compile( '\^.' )

        #
        # Open output file and write header
        #
        if f_print:
            w = open( out_file, 'w' )
            w.write( "chr\tpos\tref\tobs\tdepth\tA\ta\tC\tc\tG\tg\tT\tt\tins\t\tdel\t\tA\tC\tG\tT\tmis\ts_ratio\t0.1\tratio\t0.9\t" )
            if in_bam2:
                w.write( "obs\tdepth\tA\ta\tC\tc\tG\tg\tT\tt\tins\t\tdel\t\tA\tC\tG\tT\tmis\ts_ratio\t0.1\tratio\t0.9\tfisher\n" )

        #
        # STDIN PIPE
        #   or
        # pysam.mpileup
        #
        if input_mpileup:
            for mpileup in open( arg.input_mpileup, 'rh' ):
                data = Pileup_out( mpileup, w, threshold, mismatch_rate, min_depth )
                if data:
                    print_data( data, w, min_depth, mismatch_rate )

        elif in_bam1 and in_bam2:
            for mpileup in pysam.mpileup( '-BQ', '0', '-d', '10000000', '-f', ref_fa, in_bam1, in_bam2 ):
                data = Pileup_out( mpileup, w, threshold, mismatch_rate, min_depth )
                if data:
                    print_data( data, w, min_depth, mismatch_rate )

        elif in_bam1:
            for mpileup in pysam.mpileup( '-BQ', '0', '-d', '10000000', '-f', ref_fa, in_bam1 ):
                data = Pileup_out( mpileup, w, threshold, mismatch_rate, min_depth )
                if data:
                    print_data( data, w, min_depth, mismatch_rate )
        else:
            for mpileup in iter( sys.stdin.readline, "" ):
                data = Pileup_out( mpileup, w, threshold, mismatch_rate, min_depth )
                if data:
                    print_data( data, w, min_depth, mismatch_rate )

        if f_print:
            w.close()

    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        logging.error( ("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) ) )

############################################################
def construct_arguments():

    #
    # Arguments
    #
    parser = argparse.ArgumentParser( description = 'Fisher mutation caller' )

    parser.add_argument( '-1', '--bam1',                help = '1st bam file ( normal )',   type = str,     default = None )
    parser.add_argument( '-2', '--bam2',                help = '2nd bam file ( disease )',  type = str,     default = None )
    parser.add_argument( '-r', '--ref_fa',              help = 'Reference FASTA',           type = str,     default = None )
    parser.add_argument( '-o', '--output',              help = 'Output text file',          type = str,     default = None )
    parser.add_argument( '-q', '--base_quality',   help = 'Base quality threshold',    type = int,     default = 15 )
    parser.add_argument( '-m', '--mismatch_rate',       help = 'Mismatch rate',             type = float,   default = 0.07 )
    parser.add_argument( '-d', '--min_depth',           help = 'Mimimum depth',             type = float,   default = 9 )
    parser.add_argument( '-i', '--input_mpileup',       help = 'Input mpileupt file',       type = str,     default = None )

    #
    # Log settings
    #
    parser.add_argument( '-g', '--log_file',  help = "Log file name", type = str, default = None )
    parser.add_argument( '-l', '--log_level', help = "Logging level", type = str, default = 'DEBUG' )

    return parser
            

############################################################
def PrintHeader( myself, arg ):
    now = datetime.now()

    logging.info( '#' * 84 )
    logging.info( '# Summary' )
    logging.info( '# Generated by {my}'.format( my = myself ) )
    logging.info( '# %(y)d.%(m)d.%(d)d.%(h)d:%(m)d' % { 'y': now.year, 'm': now.month, 'd': now.day, 'h': now.hour, 'm': now.minute } )
    logging.info( '#' * 84 + '' )
    logging.info( "bam1: {0}".format( arg.bam1 ) )
    logging.info( "bam2: {0}".format( arg.bam2 ) )
    logging.info( "output: {0}".format( arg.output ) )
    logging.info( "reference fastq: {0}".format( arg.ref_fa ) )
    logging.info( "quality_threshold: {0}".format( arg.base_quality ) )
    logging.info( "mismatch_rate: {0}".format( arg.mismatch_rate ) )
    logging.info( "min_depth: {0}".format( arg.min_depth ) )
    logging.info( '-' * 84 + '' )

############################################################
#
# Main
#
def main():

    #
    # Arguments parse
    #
    argvs = sys.argv
    myself = argvs[ 0 ]
    argc = len(argvs)

    arg_parser = construct_arguments()
    
    if argc < 2:
        arg_parser.print_help()
        sys.exit(1)
                                
    global arg
    arg = arg_parser.parse_args()

    #
    # logging setup
    #
    # Level     function            value    description
    # CRITICAL  logging.critical()  50      Output only critical errors
    # ERROR     logging.error()     40      Output errors
    # WARNING   logging.warning()   30      Output warnings
    # INFO      logging.info()      20      Output information
    # DEBUG     logging.debug()     10      Output debug information
    # NOTSET                        0       Output all
    #
    level = logging.getLevelName( arg.log_level )

    if arg.log_file:
        logging.basicConfig( filename   = arg.log_file,
                             level      = level,
                             format     = '%(asctime)s %(message)s',
                             datefmt    ='%m/%d/%Y %I:%M:%S%p' )
    else:
        logging.basicConfig( level      = level,
                             format     = '%(asctime)s %(message)s',
                             datefmt    ='%m/%d/%Y %I:%M:%S%p' )
    #
    # Main function
    #
    Pileup_and_count( 
            in_bam1 = arg.bam1,
            in_bam2 = arg.bam2,
            out_file = arg.output,
            input_mpileup = arg.input_mpileup,
            ref_fa = arg.ref_fa,
            threshold = arg.base_quality,
            mismatch_rate = arg.mismatch_rate,
            min_depth = arg.min_depth,
          )

################################################################################
if __name__ == '__main__':
    main()

