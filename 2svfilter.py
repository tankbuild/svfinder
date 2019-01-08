#!/usr/bin/python
import itertools
import sys
import re
import argparse
import os
# from collections import defaultdict


"""
This script require python3 to run.
This script is to convert pindel original output file into svfilter format!
We try to combine svdetect and pindel SVs and to analyze the differences bettween male and female.

This script will also filter data by excluding up|down-stream mapped reads numbers less than 3, and reads quality less than 10.
Meanwhile structural size no less than 50bp were kept. Absolutely, you can change thses parameters by give corresponding option.

Upstream read range is between maximum and mimimum of the position of the mapped half of all upstream paired-end reads -+ half read length.
So is downstream read range.

svdetect input file should start with sex(female or male) and end with filtered.
"""


class Parser():

    def __init__(self):
        parser = argparse.ArgumentParser(
            description='Convert original output files of different SVs detection software into svfilter acceptable format files',
            usage='''python3 2svfilter.py <command> [options]

Command:  pindel\tConvert pindel output original output files into svfilter acceptable format files
\t  svdetect\tConvert svdetect output original output files into svfilter acceptable format files
'''
        )
        parser.add_argument('command', help='Command to run', nargs='?')
        args = parser.parse_args(sys.argv[1:2])
        if not args.command:
            print()
            parser.print_usage()
            print()
            exit(1)
        if not hasattr(self, args.command):
            print()
            parser.print_usage()
            print()
            exit(1)
        getattr(self, args.command)()

    def pindel(self):
        parser = argparse.ArgumentParser(
            description='Convert pindel original output file into bed and svfilter ',
            usage='''python3 radseq_analysis.py pindel -i input_folder [-q] [-n] [-l]

Options:  -i\t--input-folder\tPath to a folder containing pindel original output file
\t  -q\t--reads-q\tset a threshold to filter out reads quality
\t  -n\t--reads-n\tset a threshold for up|downstream mapped reads numbers
\t  -l\t--read-len\tset a read length
\t  -s\t--sv-size\tset a structural size for filtering
''')
        parser.add_argument('--input-folder', '-i',
                            help='Path to a folder containing pindel original output file')
        parser.add_argument('--reads-q', '-q',
                            help='set a threshold to filter out reads quality', nargs='?', const=10, default=10, type=int)
        parser.add_argument('--reads-n', '-n',
                            help='set a threshold for up|downstream mapped reads numbers', nargs='?', const=3, default=3, type=int)
        parser.add_argument('--read-len', '-l',
                            help='set a read length', nargs='?', const=150, default=150, type=int)
        parser.add_argument('--sv-size', '-s',
                            help='set a structural size for filtering', nargs='?', const=50, default=50, type=int)
        # parser.add_argument('--output-file', '-o',
        #                     help='Path to output file', nargs='?',
        #                     default='haplotypes_matrix.tsv')
        args = parser.parse_args(sys.argv[2:])
        if not args.input_folder or not os.path.isdir(args.input_folder):
            print('\nError: no valid input folder specified\n')
            parser.print_usage()
            print()
            exit(1)
        print('reads numbers: {}'.format(args.reads_n))
        print('reads quality: {}'.format(args.reads_q))
        print('read length: {}'.format(args.read_len))
        print('sv size: {}'.format(args.sv_size))
        analysis(input_folder=args.input_folder,
                 reads_n=args.reads_n,
                 reads_q=args.reads_q,
                 read_len=args.read_len,
                 sv_size=args.sv_size,
                 analysis='pindel')
        # if not args.popmap or not os.path.isfile(args.popmap):
        #     print('\nError: no valid popmap file specified\n')
        #     parser.print_usage()
        #     print()
        #     exit(1)

    def svdetect(self):
        parser = argparse.ArgumentParser(
            description='Convert svdetect original output file into bed and svfilter ',
            usage='''python3 radseq_analysis.py svdetect -i input_folder [-n]

Options:  -i\t--input-folder\tPath to a folder containing svdetect original output file
\t  -n\t--reads-n\tset a threshold to filter out reads quality
''')
        parser.add_argument('--input-folder', '-i',
                            help='Path to a folder containing pindel original output file')
        parser.add_argument('--reads-n', '-n',
                            help='set a threshold up|downstream mapped reads', nargs='?', const=3, default=3, type=int)
        # parser.add_argument('--output-file', '-o',
        #                     help='Path to output file', nargs='?',
        #                     default='haplotypes_matrix.tsv')
        args = parser.parse_args(sys.argv[2:])
        print('reads numbers: {}'.format(args.reads_n))
        if not args.input_folder or not os.path.isdir(args.input_folder):
            print('\nError: no valid input folder specified\n')
            parser.print_usage()
            print()
            exit(1)
        analysis(input_folder=args.input_folder,
                 reads_n=args.reads_n,
                 analysis='svdetect')


def isa_group_separator(line):
    return line == '####################################################################################################\n'


def pindel2svfilter(input_folder=None,
                    reads_n=None,
                    reads_q=None,
                    read_len=None,
                    sv_size=None):
    for filename in os.listdir(input_folder):
        if filename.endswith(('_D', '_INV', '_SI', '_TD')) and os.stat(os.path.join(input_folder, filename)).st_size != 0:
            print(filename)
            file = open(os.path.join(input_folder, filename), 'r')
            filtername = filename + '.filter'
            bedname = filename + '.bed'
            svfilter_file = open(os.path.join(input_folder, filtername), 'w')
            bed_file = open(os.path.join(input_folder, bedname), 'w')
            for key, group in itertools.groupby(file, isa_group_separator):
                # print(key,list(group))  # uncomment to see what itertools.groupby does.
                if not key:
                    # print('yes')
                    chrom_name = ''
                    upstream_pos = []
                    downstream_pos = []
                    reads_number = 0
                    reads_ID = []
                    svtype = ''
                    up_reads = 0
                    down_reads = 0
                    bp_start = ''
                    bp_end = ''
                    sv_s = 0
                    for item in group:
                        if re.match(r'^\d', item):
                            # print(info)
                            split = re.split(r'\s+', item)[:-1]
                            bp_start = split[9]
                            bp_end = split[10]
                            chrom_name = split[7]
                            sv_s = int(split[2])
                            # print(split[2], split[3])
                            if split[1] == 'D':
                                # print(split)
                                svtype = 'DELETION'
                            elif split[1] == 'INV':
                                # print(split)
                                svtype = 'IVERSION'
                            elif split[1] == 'I':
                                # print(split)
                                # print(sv_s)
                                svtype = 'SHORT_INSERT'
                            elif split[1] == 'TD':
                                # print(split)
                                svtype = 'TANDEM_DUPLI'
                        if svtype == 'DELETION':
                            if re.match(r'^\w', item, re.I):
                                # print('ref')
                                split = re.split(r'\s+', item)[:-1]
                                pass
                            elif re.match(r'^\s+[A-Z]', item, re.I):
                                # print('reads_info')
                                split = re.split(r'\s+', item)[1:-1]
                                # print(split[4])
                                if re.match(r'\d+', split[4]):
                                    if int(split[4]) >= reads_q:
                                        reads_number += 1
                                        if split[2] == '+':
                                            up_reads += 1
                                            reads_ID.append(split[6][1:-2])
                                            upstream_pos.extend([int(split[3]) + read_len / 2, int(split[3]) - read_len / 2])
                                        elif split[2] == '-':
                                            down_reads += 1
                                            reads_ID.append(split[6][1:-2])
                                            downstream_pos.extend([int(split[3]) + read_len / 2, int(split[3]) - read_len / 2])
                                else:
                                    if int(split[3]) >= reads_q:
                                        reads_number += 1
                                        if split[1] == '+':
                                            up_reads += 1
                                            reads_ID.append(split[5][1:-2])
                                            upstream_pos.extend([int(split[2]) + read_len / 2, int(split[2]) - read_len / 2])
                                        elif split[1] == '-':
                                            down_reads += 1
                                            reads_ID.append(split[5][1:-2])
                                            downstream_pos.extend([int(split[2]) + read_len / 2, int(split[2]) - read_len / 2])
                        if svtype == 'SHORT_INSERT':
                            if re.match(r'^\w', item, re.I):
                                # print('ref')
                                split = re.split(r'\s+', item)[:-1]
                                pass
                            elif re.match(r'^\s+[A-Z]', item, re.I):
                                # print('reads_info')
                                split = re.split(r'\s+', item)[1:-1]
                                # print(split[1], split[2], split[3], split[5])
                                if int(split[3]) >= reads_q:
                                    # print(split[1], split[2], split[3], split[5])
                                    reads_number += 1
                                    if split[1] == '+':
                                        up_reads += 1
                                        reads_ID.append(split[5][1:-2])
                                        upstream_pos.extend([int(split[2]) + read_len / 2, int(split[2]) - read_len / 2])
                                        # print(str(int(split[2]) + read_len / 2))
                                    elif split[1] == '-':
                                        down_reads += 1
                                        reads_ID.append(split[5][1:-2])
                                        downstream_pos.extend([int(split[2]) + read_len / 2, int(split[2]) - read_len / 2])
                                    # else:
                                    #     print('out of consideration')
                        if svtype == 'IVERSION':
                            if re.match(r'^\w', item, re.I):
                                # print('ref')
                                split = re.split(r'\s+', item)[:-1]
                                pass
                            elif re.match(r'^\s+[A-Z]', item, re.I):
                                # print('reads_info')
                                split = re.split(r'\s+', item)[1:-1]
                                # print(split[2], split[3], split[5])
                                if int(split[3]) >= reads_q:
                                    reads_number += 1
                                    if split[1] == '+':
                                        up_reads += 1
                                        reads_ID.append(split[5][1:-2])
                                        upstream_pos.extend([int(split[2]) + read_len / 2, int(split[2]) - read_len / 2])
                                    elif split[1] == '-':
                                        down_reads += 1
                                        reads_ID.append(split[5][1:-2])
                                        downstream_pos.extend([int(split[2]) + read_len / 2, int(split[2]) - read_len / 2])
                                    # else:
                                    #     print('out of consideration')
                        if svtype == 'TANDEM_DUPLI':
                            if re.match(r'^\w', item, re.I):
                                # print('ref')
                                split = re.split(r'\s+', item)[:-1]
                                pass
                            elif re.match(r'^\s+[A-Z]', item, re.I):
                                # print('reads')
                                split = re.split(r'\s+', item)[1:-1]
                                pass
                            elif re.match(r'^\s+[+,-]', item):
                                # print('reads_info')
                                split = re.split(r'\s+', item)[1:-1]
                                # print(split)
                                if int(split[2]) >= reads_q:
                                    reads_number += 1
                                    if split[0] == '+':
                                        up_reads += 1
                                        reads_ID.append(split[4][1:-2])
                                        # print(type(read_len))
                                        # print(type(reads_n))
                                        # print(type(reads_q))
                                        # print(type(sv_size))
                                        # print(type(int(split[1])))
                                        upstream_pos.extend([int(split[1]) + read_len / 2, int(split[1]) - read_len / 2])
                                    elif split[0] == '-':
                                        down_reads += 1
                                        reads_ID.append(split[4][1:-2])
                                        downstream_pos.extend([int(split[1]) + read_len / 2, int(split[1]) - read_len / 2])
                                    # else:
                                    #     print('out of consideration')
                            else:
                                pass
                    # print(down_reads, up_reads, sv_s, reads_q, reads_n, read_len, sv_size)
                    if down_reads >= reads_n and up_reads >= reads_n and sv_s >= sv_size:
                        # pass
                        reads_id = ','.join(reads_ID)
                        upstream_pos_s = str(min(upstream_pos))[:-2]
                        upstream_pos_e = str(max(upstream_pos))[:-2]
                        downstream_pos_s = str(min(downstream_pos))[:-2]
                        downstream_pos_e = str(max(downstream_pos))[:-2]
                        # print(reads_id, upstream_pos_s, upstream_pos_e, downstream_pos_s, downstream_pos_e)
                        svfilter_file.write(chrom_name + '\t' + upstream_pos_s + '\t' + upstream_pos_e + '\t' + 'F' + '\t' + chrom_name + '\t' + downstream_pos_s + '\t' + downstream_pos_e + '\t' + 'R' + '\t' + str(reads_number) + '\t' + reads_id + '\t' + svtype + '\n')
                        bed_file.write(chrom_name + '\t' + bp_start + '\t' + bp_end + '\t' + svtype + '\t' + str(sv_s) + '\n')
                    else:
                        pass
                        # print(down_reads, up_reads, sv_s)


def svdetect2svfilter(input_folder=None,
                      reads_n=None):
    for filename in os.listdir(input_folder):
        if filename.endswith(".filtered"):
            if filename.startswith("male"):
                suffix = "male"
            elif filename.startswith("female"):
                suffix = "female"
            deletion_file = open("svdetect_deletion." + suffix, 'w')
            large_dupli_file = open("svdetect_large_dupli." + suffix, 'w')
            inversion_file = open("svdetect_inversion." + suffix, 'w')
            inv_dupli_file = open("svdetect_inv_dupli." + suffix, 'w')
            duplication_file = open("svdetect_duplication." + suffix, 'w')
            transloc_file = open("svdetect_transloc." + suffix, 'w')
            inv_transloc_file = open("svdetect_inv_transloc." + suffix, 'w')
            small_dupli_file = open("svdetect_small_dupli." + suffix, 'w')
            file = open(input_folder + filename, 'r')
            for line in file:
                info = re.split(r'\s+', line[:-1])
                # print(info)
                if re.match(r'(\d).*?', info[18]) and int(re.match(r'(\d).*?', info[18]).group(1)) >= reads_n:
                    if info[16] == "DELETION":
                        if re.match(r'\(R.*?', info[8]):
                            deletion_file.write(info[0] + '\t' + info[1] + '\t' + info[2] + '\t' + 'R' + '\t' + info[3] + '\t' + info[4] + '\t' + info[5] + '\t' + 'F' + '\t' + info[6] + '\t' + info[7][1:-1] + '\t' + info[16] + '\n')
                        else:
                            deletion_file.write(info[0] + '\t' + info[1] + '\t' + info[2] + '\t' + 'F' + '\t' + info[3] + '\t' + info[4] + '\t' + info[5] + '\t' + 'R' + '\t' + info[6] + '\t' + info[7][1:-1] + '\t' + info[16] + '\n')
                    elif info[16] == "INVERSION":
                        if re.match(r'\(R.*?', info[8]):
                            inversion_file.write(info[0] + '\t' + info[1] + '\t' + info[2] + '\t' + 'R' + '\t' + info[3] + '\t' + info[4] + '\t' + info[5] + '\t' + 'R' + '\t' + info[6] + '\t' + info[7][1:-1] + '\t' + info[16] + '\n')
                        else:
                            inversion_file.write(info[0] + '\t' + info[1] + '\t' + info[2] + '\t' + 'F' + '\t' + info[3] + '\t' + info[4] + '\t' + info[5] + '\t' + 'F' + '\t' + info[6] + '\t' + info[7][1:-1] + '\t' + info[16] + '\n')
                    elif info[16] == "LARGE_DUPLI":
                        if re.match(r'\(R.*?', info[8]):
                            large_dupli_file.write(info[0] + '\t' + info[1] + '\t' + info[2] + '\t' + 'R' + '\t' + info[3] + '\t' + info[4] + '\t' + info[5] + '\t' + 'F' + '\t' + info[6] + '\t' + info[7][1:-1] + '\t' + info[16] + '\n')
                        else:
                            large_dupli_file.write(info[0] + '\t' + info[1] + '\t' + info[2] + '\t' + 'F' + '\t' + info[3] + '\t' + info[4] + '\t' + info[5] + '\t' + 'R' + '\t' + info[6] + '\t' + info[7][1:-1] + '\t' + info[16] + '\n')
                    elif info[16] == "INV_DUPLI":
                        if re.match(r'\(R.*?', info[8]):
                            inv_dupli_file.write(info[0] + '\t' + info[1] + '\t' + info[2] + '\t' + 'R' + '\t' + info[3] + '\t' + info[4] + '\t' + info[5] + '\t' + 'R' + '\t' + info[6] + '\t' + info[7][1:-1] + '\t' + info[16] + '\n')
                        else:
                            inv_dupli_file.write(info[0] + '\t' + info[1] + '\t' + info[2] + '\t' + 'F' + '\t' + info[3] + '\t' + info[4] + '\t' + info[5] + '\t' + 'F' + '\t' + info[6] + '\t' + info[7][1:-1] + '\t' + info[16] + '\n')
                    elif info[16] == "DUPLICATION":
                        if re.match(r'\(R.*?', info[8]):
                            duplication_file.write(info[0] + '\t' + info[1] + '\t' + info[2] + '\t' + 'R' + '\t' + info[3] + '\t' + info[4] + '\t' + info[5] + '\t' + 'R' + '\t' + info[6] + '\t' + info[7][1:-1] + '\t' + info[16] + '\n')
                        else:
                            duplication_file.write(info[0] + '\t' + info[1] + '\t' + info[2] + '\t' + 'F' + '\t' + info[3] + '\t' + info[4] + '\t' + info[5] + '\t' + 'F' + '\t' + info[6] + '\t' + info[7][1:-1] + '\t' + info[16] + '\n')

                    elif info[16] == "TRANSLOC":
                        if re.match(r'\(R.*?', info[8]):
                            transloc_file.write(info[0] + '\t' + info[1] + '\t' + info[2] + '\t' + 'R' + '\t' + info[3] + '\t' + info[4] + '\t' + info[5] + '\t' + 'F' + '\t' + info[6] + '\t' + info[7][1:-1] + '\t' + info[16] + '\n')
                        else:
                            transloc_file.write(info[0] + '\t' + info[1] + '\t' + info[2] + '\t' + 'F' + '\t' + info[3] + '\t' + info[4] + '\t' + info[5] + '\t' + 'R' + '\t' + info[6] + '\t' + info[7][1:-1] + '\t' + info[16] + '\n')
                    elif info[16] == "SMALL_DUPLI":
                        if re.match(r'\(R.*?', info[8]):
                            small_dupli_file.write(info[0] + '\t' + info[1] + '\t' + info[2] + '\t' + 'R' + '\t' + info[3] + '\t' + info[4] + '\t' + info[5] + '\t' + 'F' + '\t' + info[6] + '\t' + info[7][1:-1] + '\t' + info[16] + '\n')
                        else:
                            small_dupli_file.write(info[0] + '\t' + info[1] + '\t' + info[2] + '\t' + 'F' + '\t' + info[3] + '\t' + info[4] + '\t' + info[5] + '\t' + 'R' + '\t' + info[6] + '\t' + info[7][1:-1] + '\t' + info[16] + '\n')
                elif re.match(r'(\d).*?', info[17]) and int(re.match(r'(\d).*?', info[17]).group(1)) >= reads_n:
                    if info[16] == "INV_TRANSLOC":
                        if re.match(r'\(R.*?', info[8]):
                            inv_transloc_file.write(info[0] + '\t' + info[1] + '\t' + info[2] + '\t' + 'R' + '\t' + info[3] + '\t' + info[4] + '\t' + info[5] + '\t' + 'R' + '\t' + info[6] + '\t' + info[7][1:-1] + '\t' + info[16] + '\n')
                        else:
                            inv_transloc_file.write(info[0] + '\t' + info[1] + '\t' + info[2] + '\t' + 'F' + '\t' + info[3] + '\t' + info[4] + '\t' + info[5] + '\t' + 'F' + '\t' + info[6] + '\t' + info[7][1:-1] + '\t' + info[16] + '\n')


def analysis(input_folder=None,
             reads_n=None,
             reads_q=None,
             read_len=None,
             sv_size=None,
             analysis=None):
    if analysis == 'pindel':
        pindel2svfilter(input_folder, reads_n, reads_q, read_len, sv_size)
    elif analysis == 'svdetect':
        svdetect2svfilter(input_folder, reads_n)


Parser()
