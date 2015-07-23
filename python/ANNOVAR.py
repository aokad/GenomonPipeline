#! /usr/bin/python
"""

This is a Python parser module for Annovar output file


"""
import sys
import os
import csv

#
# Header
#
##FILTER=<ID=LowQual,Description="QUAL < 50.0">
##INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP Membership">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=DS,Number=0,Type=Flag,Description="Were any of the samples downsampled?">
##INFO=<ID=Dels,Number=1,Type=Float,Description="Fraction of Reads Containing Spanning Deletions">
##INFO=<ID=HRun,Number=1,Type=Integer,Description="Largest Contiguous Homopolymer Run of Variant Allele In Either Direction">
##INFO=<ID=HaplotypeScore,Number=1,Type=Float,Description="Consistency of the site with two (and only two) segregating haplotypes">
##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">
##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
##INFO=<ID=SB,Number=1,Type=Float,Description="Strand Bias">
##INFO=<ID=VQSLOD,Number=1,Type=Float,Description="log10-scaled probability of variant being true under the trained gaussian mixture model">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth (only filtered reads used for calling)">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=3,Type=Float,Description="Normalized, Phred-scaled likelihoods for AA,AB,BB genotypes where A=ref and B=alt; not applicable if site is not biallelic">
#CHROM  POS ID      REF ALT QUAL    FILTER  INFO    FORMAT
#
# INFO example
# AC=1;AF=0.50;AN=2;DP=315;Dels=0.00;HRun=2;HaplotypeScore=15.11;MQ=91.05;MQ0=15;QD=16.61;SB=-1533.02;VQSLOD=-1.5473
#
# FORMAT example
# GT:AD:DP:GQ:PL   0/1:173,141:282:99:255,0,255
#
# Annotation tag in VCF     Meaning
# AC,AF,AN    See the Technical Documentation for Chromosome Counts.
# DB  If present, then the variant is in dbSNP.
# DP  See the Technical Documentation for Coverage.
# DS  Were any of the samples downsampled because of too much coverage?
# Dels    See the Technical Documentation for SpanningDeletions.
# MQ and MQ0  See the Technical Documentation for RMS Mapping Quality and Mapping Quality Zero.
# BaseQualityRankSumTest  See the Technical Documentation for Base Quality Rank Sum Test.
# MappingQualityRankSumTest   See the Technical Documentation for Mapping Quality Rank Sum Test.
# ReadPosRankSumTest  See the Technical Documentation for Read Position Rank Sum Test.
# HRun    See the Technical Documentation for Homopolymer Run.
# HaplotypeScore  See the Technical Documentation for Haplotype Score.
# QD  See the Technical Documentation for Qual By Depth.
# VQSLOD  Only present when using Variant quality score recalibration. Log odds ratio of being a true variant versus being false under the trained gaussian mixture model.
# FS  See the Technical Documentation for Fisher Strand
# SB  How much evidence is there for Strand Bias (the variation being seen on only the forward or only the reverse strand) in the reads? Higher SB values denote more bias (and therefore are more likely to indicate false positive calls).

################################################################################
#
# ANNOVAR output parser
#
################################################################################
class AnnoCSVOut:

    def __init__( self, filename ):
        self.data = {}

        self.data[ 'filename' ] = filename;

        try:
            self.file = open( filename, 'r' )

        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print( ("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) ) )
            raise

    def __del__( self ):
        try:
            if self.file:
                self.file.close()

        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print( ("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) ) )
            raise

    def close( self ):
        self.__del__()

    def read( self ):

        #
        # Phred score 0 ~ 41
        qual = 41

        #
        # PASS,
        filter = 'PASS'
        info = '.'
        format = '.'

        try:
            csv_obj = csv.reader( self.file )
            header = csv_obj.next();
            id = 0
            header_id = {}
            for item in header:
                header_id[ item ] = id
                id += 1
                if 'dbSNP' in item:
                    dbSNP_version = item

            for line in csv_obj:
                id = line[ header_id[ dbSNP_version ] ]
                id = id if id != '' else '.'

                yield_value = (
                        line[ header_id[ 'Chr' ] ],
                        line[ header_id[ 'Start' ] ],
                        line[ header_id[ 'End' ] ],
                        id,
                        line[ header_id[ 'Ref' ] ],
                        line[ header_id[ 'Obs' ] ],
                        qual,
                        filter,
                        info,
                        format,
                        )

                yield yield_value


        except ValueError:
            print "Error parse file"

class AnnoTSVOut:

    def __init__( self, filename ):
        self.data = {}

        self.data[ 'filename' ] = filename;

        try:
            self.file = open( filename, 'r' )

        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print( ("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) ) )
            raise

    def __del__( self ):
        try:
            if self.file:
                self.file.close()

        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print( ("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) ) )
            raise

    def close( self ):
        self.__del__()

    def read( self ):

        #
        # Phred score 0 ~ 41
        qual = 41

        #
        # PASS,
        filter = 'PASS'
        info = '.'
        format = '.'

        try:
            header = self.file.readline().split( "\t" )
            id = 0
            header_id = {}
            for item in header:
                header_id[ item ] = id
                id += 1
                if 'dbSNP' in item:
                    dbSNP_version = item

            for line_tmp in self.file:
                line = line_tmp.split( "\t" )
                id = line[ header_id[ dbSNP_version ] ]
                id = id if id != '' else '.'

                yield_value = (
                        line[ header_id[ 'Chr' ] ],
                        line[ header_id[ 'Start' ] ],
                        line[ header_id[ 'End' ] ],
                        id,
                        line[ header_id[ 'Ref' ] ],
                        line[ header_id[ 'Obs' ] ],
                        qual,
                        filter,
                        info,
                        format
                        )

                yield yield_value


        except ValueError:
            print "Error parse file"

################################################################################
#
# ANNOVAR input file parser
#
################################################################################
class AnnoIn:

    def __init__( self, filename ):
        self.data = {}

        self.data[ 'filename' ] = filename;

        try:
            self.file = open( filename, 'r' )

        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print( ("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) ) )
            raise

    def __del__( self ):
        try:
            if self.file:
                self.file.close()

        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print( ("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) ) )
            raise

    def close( self ):
        self.__del__()

    def read( self ):

        #
        # Phred score 0 ~ 41
        qual = 41

        #
        # PASS,
        filter = 'PASS'
        info = '.'
        format = '.'

        try:

            #
            # ANNOVAR output format
            # <barcode> chr start end ref obs depth A,C,G,T mis_rate strand_ratio 0.1 ratio 0.9
            # <fisher>  chr start end ref obs depth A,C,G,T tumor_mis_rate tumor_strand_ratio control_mis_rate control_strand_ratio fisher
            #
            POS_CHR     = 0
            POS_START   = 1
            POS_END     = 2
            POS_REF     = 3
            POS_OBS     = 4
            POS_DEPTH   = 5
            for line in self.file:
                if line[ 0 ] ==  '#':
                    continue

                line_split = line.split( '\t' )
                yield_value = (
                        line_split[ POS_CHR ],
                        line_split[ POS_START ],
                        line_split[ POS_END ],
                        '.',
                        line_split[ POS_REF ],
                        line_split[ POS_OBS ],
                        qual,
                        filter,
                        info,
                        format
                        )

                yield yield_value

        except ValueError:
            print "Error parse file"

