#!/usr/bin/env python
# usage: get_enriched_regions.py [OPTIONS] bam_file
# v0.1a
import os
import sys
import argparse
import collections
import statistics
import os.path as path
import pysam


SELFPATH = path.dirname(__file__)
parser_ar = argparse.ArgumentParser(prog='get_enriched_regions',
                                    usage='get_enriched_regions [OPTION] <bam_file>',
                                    description='Get enriched regions for a bam file',
                                    epilog='',
                                    formatter_class=argparse.RawTextHelpFormatter)

parser_ar.add_argument('bam_file', help='FILE. A sorted and indexed bam file', metavar='BAM')

parser_ar.add_argument('-o', help='FOLDER. The Folder to store the statistics results [bam file folder]', metavar='[folder]', dest='OUTPUT')
parser_ar.add_argument('-m', metavar='[50]', default=50, type=int, help='INT. Minimum coverage for the enriched regions.', dest='MINIMUM_COVERAGE')
parser_ar.add_argument('-n', metavar='[1]', default=3, type=int, help='INT. Enriched rigions must have more than n base with minimum coverage.', dest='MINIMUM_N_BASE')
parser_ar.add_argument('-g', metavar='[3]', default=3, type=int, help='INT. Minimum gap between two adjacent enriched regions.', dest='MAXIMUM_MAP_REGION_GAP')
parser_ar.add_argument('-v', action='store_true', default=False, help="Print out the error messages if there's any.", dest='IS_VERBOSE')


paramters = parser_ar.parse_args()

BAM_FILE = paramters.bam_file

OUTPUT = paramters.OUTPUT
MINIMUM_COVERAGE = paramters.MINIMUM_COVERAGE
N_BASES = paramters.MINIMUM_N_BASE
MAXIMUM_MAP_REGION_GAP = paramters.MAXIMUM_MAP_REGION_GAP
IS_VERBOSE = paramters.IS_VERBOSE

# Based on some customed rules, get qualified regions,


def get_region(bamfile: str, minimum_coverage: int, n_bases: int, maximum_gap_length: int, verbose = False) -> list[tuple]:
    '''
    This function will return all qualified regions according to minimum_coverage, maximum_gap_length
    The region must have at least n_bases with higher (or equal) coverage than minimum_coverage and lower gap length than maximum_gap_length

    Parameters:
        **bamfile**: str
            The bam file to be calculated. We assume it has been coordinate-sorted and indexed.

        **minimum_coverage**: type
            its function

        **n_bases**: type
            its function

        **maximum_gap_length**: type
            its function

    Returns:
        **value**: list[tuple]
            its meaning
    '''
    region_list = []  # store the result

    is_qualified = False
    current_chrom_str = ''
    current_int = -1

    start_int = -1
    end_int = -1
    end_chrom_str = ''

    n = 0
    bamfile_af = pysam.AlignmentFile(bamfile, 'rb')
    for pileupcolumn in bamfile_af.pileup(stepper = 'all', min_base_quality= 0, ignore_overlaps= False, ignore_orphans= False):
        if current_chrom_str != pileupcolumn.reference_name:
            print('Searching enriched regions at contig {}...'.format(pileupcolumn.reference_name))

        current_chrom_str = pileupcolumn.reference_name
        current_int = pileupcolumn.reference_pos
        if start_int == -1:
            start_int = current_int
            end_int = current_int
            end_chrom_str = current_chrom_str

        try:
            if pileupcolumn.nsegments >= minimum_coverage:
                n += 1
                if n >= n_bases:
                    is_qualified = True
        except:
            message = 'Fail to get coverage number at "{chr}:{pos}". Skip'.format(chr= pileupcolumn.reference_name, pos= pileupcolumn.reference_pos)
            if IS_OUTPUT_ERROR:
                print(message)
            continue

        # on following situation we form a region
        region_tu = ()
        if current_int - end_int > maximum_gap_length + 1 or end_chrom_str != pileupcolumn.reference_name:
            region_tu = (end_chrom_str, start_int, end_int)
            start_int = current_int
            end_int = current_int
            end_chrom_str = current_chrom_str
        else:
            end_int = current_int

        # should we add
        if is_qualified and region_tu != ():
            region_list.append(region_tu)
            if verbose:
                print(region_tu)
            is_qualified = False
            n = 0

    bamfile_af.close()

    return region_list

# Given a region, extract some statisic (eg. ['width', 'mean_coverage', 'mean_MAPQ', 'total_reads', 'specific_primers', 'mean_aligned_length'])


def extract_info(bamfile: str, region: tuple, primer_length = 20) -> dict:
    '''
    This function will extract some statisic given a region
    The specific primer sequenced are extracted for read1 5'end based on *primer_length*

    Parameters:
        **bamfile**: str
            The bam file to be calculated. We assume it has been coordinate-sorted and indexed.

        **region**: tuple
            its function

        **primer_length**: type
            its function

    Returns:
        **value**: dict
            its meaning
    '''

    result_dict = {}
    result_dict['contig'] = region[0]
    result_dict['start'] = region[1]
    result_dict['end'] = region[2]
    result_dict['width'] = region[2] - region[1] + 1

    coverage_list = []
    bamfile_af = pysam.AlignmentFile(bamfile, 'rb')
    for pileupcolumn in bamfile_af.pileup(contig = region[0], start = region[1], stop = region[2], stepper = 'all', min_base_quality= 0, ignore_overlaps= False, ignore_orphans= False):
        coverage_list.append(pileupcolumn.nsegments)

    result_dict['mean_coverage'] = round(statistics.mean(coverage_list))

    map_qualtity_list = []
    specific_primers = []  # at read1 5'-end
    universial_primers = []  # universial primer is not sequenced so this is only for counting read2
    aligned_length = []
    for segment in bamfile_af.fetch(contig = region[0], start = region[1], stop = region[2]):
        map_qualtity_list.append(segment.mapping_quality)
        aligned_length.append(len(segment.query_alignment_sequence))
        primer_str = segment.query_alignment_sequence[:primer_length]
        if segment.is_read1:
            specific_primers.append(primer_str)
        elif segment.is_read2:
            universial_primers.append(primer_str)
        else:
            message = 'Read {name} at {chrom}:{start}-{end} is neither read1 nor read2. Skip'.format(segment.query_name, segment.reference_name, segment.reference_start, segment.reference_end)
            print(message)
            continue

    result_dict['mean_MAPQ'] = round(statistics.mean(map_qualtity_list))
    result_dict['total_reads'] = len(specific_primers) + len(universial_primers)
    result_dict['mean_aligned_length'] = round(statistics.mean(aligned_length))
    result_dict['specific_primers'] = collections.Counter(specific_primers)

    bamfile_af.close()
    return result_dict


def output_info(result: list[dict], file = 'enriched_region.csv') -> int:

    with open(file, 'wt') as out_f:
        line_str = 'chromosome\tstart\tend\twidth\tcoverage\tMAPQ\tmean_aligned_length\ttotal_reads\tspecific_primers\n'
        out_f.writelines(line_str)
        for result_dict in result:
            '''
            {'contig': 'chr19',
            'start': 54552521,
            'end': 54552853,
            'width': 333,
            'mean_coverage': 30,
            'mean_MAPQ': 58,
            'total_reads': 118,
            'mean_aligned_length': 80,
            'specific_primers': Counter({'GTGTGTGTATGTGAGTGTGT': 20,
                              'GTGTGTGCATGTGAGATTGT': 17,
                              'TGTGTGTATGTGAGTGTGTG': 2,
                              'TGTGAGTGTGTGAGATTGTA': 2,
                              'TGCATGTGAGATTGTATGAG': 2,
                              'TGTGTGTGTGTGAGTTTGAG': 1,
                              'GTGTGAGTATATGTATGTGT': 1,
                              'TGTGTATGTGTGTGTATGTG': 1,
                              'ATGTGTGTATGTGAGTGTGT': 1,
                              'TGTATGTGAGTGTGTGAGAT': 1,
                              'TGAGATTGTATGTGTGTGCA': 1,
                              'GTGTGCATGTTGTGTGTGTC': 1,
                              'TGTTGTGTGTGTCTGAGATT': 1,
                              'TTGTGTGTGTGTGTGTGCAT': 1,
                              'TGTGTGTGTGTGTGTGCATG': 1,
                              'TGTGTGCATGTGAGATTGTA': 1,
                              'GTGCATGTGAGATTGTATGA': 1,
                              'ATGTGAGATTGTATGAGTGT': 1,
                              'TGTGATATTGTATGAGTGTT': 1,
                              'TGAGATTGTATGAGTGTTTG': 1,
                              'TATGAGTGTTTGCATGGTGT': 1})}
            '''

            temp_dict = result_dict['specific_primers']
            temp_list = temp_dict.most_common(len(temp_dict))  # a list of tuples
            specific_primers_str = ''
            for tu in temp_list:
                specific_primers_str += '{seq}({count}), '.format(seq = tu[0], count = tu[1])

            specific_primers_str = specific_primers_str.strip(', ')

            line_str = '{contig}\t{start}\t{end}\t{width}\t{coverage}\t{MAPQ}\t{mean_aligned_length}\t{total_reads}\t{specific_primers}\n'.format(contig = result_dict['contig'],
                                                                                                                                                  start = result_dict['start'],
                                                                                                                                                  end = result_dict['end'],
                                                                                                                                                  width = result_dict['width'],
                                                                                                                                                  coverage = result_dict['mean_coverage'],
                                                                                                                                                  MAPQ = result_dict['mean_MAPQ'],
                                                                                                                                                  total_reads = result_dict['total_reads'],
                                                                                                                                                  mean_aligned_length = result_dict['mean_aligned_length'],
                                                                                                                                                  specific_primers = specific_primers_str)
            out_f.writelines(line_str)

    return 0


def main(argvList=sys.argv, argv_int=len(sys.argv)):
    '''
    This is the main function

    Returns:
        **value**: type
            its meaning
    '''
    bam_file = BAM_FILE
    bam_file = path.realpath(path.expanduser(bam_file))
    out_file = OUTPUT
    if out_file is None:
        out_file = path.splitext(bam_file)[0] + '_regions.csv'

    region_list = get_region(bam_file,
                             minimum_coverage = MINIMUM_COVERAGE,
                             n_bases = N_BASES,
                             maximum_gap_length = MAXIMUM_MAP_REGION_GAP,
                             verbose = IS_VERBOSE)

    result_list = []
    i = 0
    for region in region_list:
        i += 1

        print()
        result_dict = extract_info(bam_file, region)
        result_list.append(result_dict)

    r = output_info(result_list, file = out_file)
    print('Done({})'.format(r))

    return


main()

