import sys
import os
import re
import pysam
import scipy.special
import fisher

class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
        return value


def Pileup_and_count( in_bam1, in_bam2, out_file, ref_fa, threshold, mismatch_rate ):

    f_print = True

    try:
        if f_print:
            w = open( out_file, 'w' )

        filter_quals = ''
        for qual in range( 33, 33 + threshold ):
            filter_quals += str( unichr( qual ) )

        if f_print:
            w.write( "chr\tpos\tref\tobs\tdepth\tA\ta\tC\tc\tG\tg\tT\tt\tins\t\tdel\t\tA\tC\tG\tT\tmis\ts_ratio\t0.1\tratio\t0.9\t" )
            w.write( "obs\tdepth\tA\ta\tC\tc\tG\tg\tT\tt\tins\t\tdel\t\tA\tC\tG\tT\tmis\ts_ratio\t0.1\tratio\t0.9\n" )

        #
        # Setup regular expression
        # ([\+\-])[0-9]+[ACGTNacgtn]+
        #
        target = re.compile( '([\+\-])([0-9]+)([ACGTNacgtn]+)' )
        remove_chr = re.compile( '\^.' )

        if in_bam2:
            mpileups = pysam.mpileup( '-BQ', '0', '-d', '10000000', '-f', ref_fa, in_bam1, in_bam2 )
        else:
            mpileups = pysam.mpileup( '-BQ', '0', '-d', '10000000', '-f', ref_fa, in_bam1)

        for mpileup in mpileups:
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
            mp_list = str( mpileup.translate( None, '\n' ) ).split( '\t' )
            ref_base = mp_list[ 2 ].upper()
            coordinate = mp_list[ 0:3 ]

            if mp_list[ 3 ] == 0 or ( in_bam2 and mp_list[ 6 ] == 0 ): # skip if depth is 0
                continue

            if ref_base == 'N': # skip if reference is 'N'
                continue

            indel = AutoVivification()

            mis_rate_list = []
            s_rate_list = []
            outstr = ''

            if in_bam2:
                input_list = ( ( mp_list[ 3 ], mp_list[ 4 ], mp_list[ 5 ] ),
                               ( mp_list[ 6 ], mp_list[ 7 ], mp_list[ 8 ] ) )
            else:
                input_list = ( ( mp_list[ 3 ], mp_list[ 4 ], mp_list[ 5 ] ) )

            for depth, read_bases, qual_list in input_list:

                if outstr != '':
                    outstr += '\t'

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

                    key = '\t'.join( coordinate + [ read_bases ] ) + "\t" + bases.upper()
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

                if len( read_bases ) != len( qual_list ):
                    print "mpileup data is not good: {0}, {1}".format( mpileup, read_bases )
                    raise

                u_ref = ref_base
                l_ref = ref_base.lower()

                read_bases = read_bases.replace( '.', u_ref )
                read_bases = read_bases.replace( ',', l_ref )

                base2qual = {
                    "A": 0,
                    "C": 0,
                    "G": 0,
                    "T": 0,
                    "a": 0,
                    "c": 0,
                    "g": 0,
                    "t": 0,
                    "N": 0,
                    "n": 0,
                    "*": 0
                }
                for nuc, qual in zip( read_bases, qual_list ):
                    if not ( qual in filter_quals ):
                        base2qual[ nuc ] += 1

                tmp_list = [ read_bases, depth ]
                for nuc in ( 'A', 'a', 'C', 'c', 'G', 'g', 'T', 't' ):
                    tmp_list.append( str( base2qual[ nuc ] ) )
                outstr += '\t'.join( tmp_list )

                for type in ( '+', '-' ):
                    if type in indel:
                        for key in indel[ type ].keys():
                            bases = key.split( '\t' )[ 2 ]
                            outstr += "\t{bases}:{p}\t{bases}:{n}".format( 
                                           bases = bases,
                                           p = indel[ type ][ key ][ '1' ],
                                           n = indel[ type ][ key ][ '-1' ] )
                    else:
                        outstr += "\t\t"

                if depth > 9:
                    #
                    # Variation count
                    #
                    count = {}
                    mis_num = 0
                    mis_base = 'N'
                    for nuc in ( 'A', 'C', 'G', 'T' ):
                        base_num = base2qual[ nuc ] + base2qual[ nuc.lower() ]
                        count.update( { nuc: base_num } )
                        if mis_num < base_num:
                            mis_num = base_num
                            mis_base = nuc

                    #
                    # Reference count
                    #
                    ref_num = count[ ref_base ]

                    #
                    # Calculate ratio
                    #
                    mis_rate = float( mis_num / int( depth ) )

                    if mis_base != 'N' and mis_rate > mismatch_rate:
                        s_ratio = float( base2qual[ mis_base ] / count[ mis_base ] )

                        outstr += "\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(
                                    count[ 'A' ],
                                    count[ 'C' ],
                                    count[ 'G' ],
                                    count[ 'T' ],
                                    mis_rate,
                                    s_ratio )

                        # p-th quantile, mis bas count, ref base count
                        outstr += "\t{0}\t{1}\t{2}".format( 
                                    scipy.special.btdtri( mis_num + 1, ref_num + 1, 0.1 ),
                                    float( ( mis_num + 1 ) / (ref_num + mis_num + 2 ) ),
                                    scipy.special.btdtri( mis_num + 1, ref_num + 1, 0.9 ) ) 

                    else:
                        outstr += '\t\t\t\t\t\t\t\t\t'

                else:
                    outstr += '\t\t\t\t\t\t\t\t\t'

            if f_print:
                w.write( '\t'.join( coordinate ) + '\t' + outstr + '\n' )

        if f_print:
            w.close()

    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        logging.error("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) )


if __name__ == '__main__':

    #
    # 0) myself
    # 1) input bam 1
    # 2) input bam 2
    # 3) output file
    # 4) reference fasta
    # 5) quality threshold
    # 6) mismatch rate
    #

    Pileup_and_count( sys.argv[ 1 ],
            sys.argv[ 2 ],
            sys.argv[ 3 ],
            sys.argv[ 4 ],
            int( sys.argv[ 5 ] ),
            float( sys.argv[ 6 ] )
          )

