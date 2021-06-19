#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Fri June 19 2020
# @author: Chuanyi Zhang

"""
Get frequency of junction-spanning reads
"""

import pysam
import pandas as pd
import sys
import argparse
from os import listdir
from os.path import isfile, join
from collections import defaultdict, Counter
from typing import List, Tuple, Dict


def read_pair_generator(bam, region_string=None, pos=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = dict()
    for read in bam.fetch(region=region_string):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            read_dict[qname] = read
        else:
            if pos == None:
                if read.is_read1:
                    yield read, read_dict[qname]
                else:
                    yield read_dict[qname], read
            elif read.get_overlap(pos, pos+1) == 1:
                yield read, read_dict[qname]
            elif read_dict[qname].get_overlap(pos, pos+1) == 1:
                yield read_dict[qname], read
            del read_dict[qname]


class IncrementalAverage():
    def __init__(self):
        self.n = 0
        self.mean = 0

    def update(self, num):
        self.mean = (self.mean * self.n + num) / (self.n + 1)
        self.n += 1


class SpliceJunction():
    def __init__(self, reference: str, min_base_qual=20, min_mapping_qual=20, min_snv_qual=20, contig=None):
        self.CIGAR = ["M", "I", "D", "N", "S", "H", "P", "=", "X", "B"]
        self.sgRNAs = {  # 0-based positions
            "S": 21562,
            "3a": 25392,
            "E": 26244,
            "M": 26522,
            "6": 27201,
            "7a": 27393,
            "7b": 27755,
            "8": 27893,
            "N": 28273,
            "10": 29557
        }
        self.ORF1ab_start = (265, 268)
        self.LeaderStart = 54
        self.LeaderEnd = 85
        self.AfterStartCodon = 60
        if not contig:
            self.contig = "NC_045512.2"
        else:
            self.contig = contig
        self.min_base_quality = min_base_qual
        self.min_mapping_quality = min_mapping_qual
        self.min_snv_qual = min_snv_qual
        self.reference = pysam.FastaFile(reference)
        self.refLength = self.reference.lengths[0]
        self.contig = self.reference.references[0]

    def identify_ORF(self, pos) -> str:

        for key, value in self.sgRNAs.items():
            if pos == value:
                return key

        return 'non-canonical'

    @staticmethod
    def convert_query_index_to_ref(query_index, read) -> int:
        ref_start = read.reference_start
        query_start = 0
        for st, length in read.cigartuples:
            if st == 0 or st == 7 or st == 8:
                if query_start <= query_index < query_start + length:
                    return ref_start + query_index - query_start
                ref_start += length
                query_start += length
            elif st == 2 or st == 3:
                ref_start += length
            elif st == 1 or st == 4:
                if query_start <= query_index < query_start + length:
                    return None
                query_start += length
        return -1

    @staticmethod
    def convert_ref_index_to_query(ref_idx, read) -> int:
        ref_start = read.reference_start
        query_start = 0
        for st, length in read.cigartuples:
            if st == 0 or st == 7 or st == 8:
                if ref_start <= ref_idx < ref_start + length:
                    return query_start + ref_idx - ref_start
                ref_start += length
                query_start += length
            elif st == 2 or st == 3:
                if ref_start <= ref_idx < ref_start + length:
                    return None
                ref_start += length
            elif st == 1 or st == 4:
                query_start += length
        return -1

    def get_spliced_junctions(self, read) -> Tuple[List[Tuple[int, int]], List[Tuple[int, int]]]:
        sj_ref = []
        sj_query = []
        ref_start = read.reference_start
        query_start = 0
        for st, length in read.cigartuples:
            if st == 3:
                # All flanking 5 nucleotides have good base quality
                if min(read.query_qualities[max(query_start-6, 0):query_start]) >= self.min_base_quality and min(read.query_qualities[query_start:min(query_start+5, read.query_length)]) >= self.min_base_quality:
                    sj_ref.append((ref_start-1, ref_start+length))
                    sj_query.append((query_start-1, query_start))
                ref_start += length
            elif st == 0 or st == 7 or st == 8:
                ref_start += length
                query_start += length
            elif st == 2:
                ref_start += length
            elif st == 1 or st == 4:
                query_start += length
        return sj_ref, sj_query

    def get_spliced_junctions_muts(self, read, mutList=[]):
        sj_ref = []
        ref_start = read.reference_start
        query_start = 0
        for st, length in read.cigartuples:
            if st == 3:
                # All flanking 5 nucleotides have good base quality
                if read.query_qualities:
                    if (min(read.query_qualities[max(query_start-6, 0):query_start]) >= self.min_base_quality and
                            min(read.query_qualities[query_start:min(query_start+5, read.query_length)])) >= self.min_base_quality:
                        sj_ref.append(
                            (ref_start-1, ref_start+length, 'splice'))
                else:
                    sj_ref.append((ref_start-1, ref_start + length, 'splice'))
                ref_start += length
            elif st == 0 or st == 7 or st == 8:
                for pos in mutList:
                    if pos < ref_start + length and ref_start <= pos:
                        shift = pos - ref_start
                        sj_ref.append(
                            (pos, pos, read.query_sequence[query_start + shift - 1: query_start + shift]))
                ref_start += length
                query_start += length
            elif st == 2:
                ref_start += length
            elif st == 1 or st == 4:
                query_start += length
        return sj_ref, ref_start

    def get_read_segments(self, read):
        read_segments = []
        sj_ref = []
        ref_start = read.reference_start
        query_start = 0
        for st, length in read.cigartuples:
            if st == 3:
                # All flanking 5 nucleotides have good base quality
                if read.query_qualities:
                    if (min(read.query_qualities[max(query_start-6, 0):query_start]) >= self.min_base_quality and
                            min(read.query_qualities[query_start:min(query_start+5, read.query_length)])) >= self.min_base_quality:
                        sj_ref.append(
                            (ref_start-1, ref_start+length, 'splice'))
                else:
                    sj_ref.append((ref_start-1, ref_start + length, 'splice'))
                ref_start += length
            elif st == 0 or st == 7 or st == 8:
                read_segments.append((ref_start, ref_start + length - 1))
                ref_start += length
                query_start += length
            elif st == 2:
                ref_start += length
            elif st == 1 or st == 4:
                query_start += length
        return read_segments, sj_ref, ref_start

    @staticmethod
    def get_read_segments_junctions(read):
        read_segments = []
        sj_ref = []
        ref_start = read.reference_start
        query_start = 0
        for st, length in read.cigartuples:
            if st == 3:
                # skipped region from reference
                sj_ref.append((ref_start-1, ref_start + length, 'splice'))
                ref_start += length
            elif st == 0 or st == 7 or st == 8:
                # alignment match, sequence match or sequence mismatch
                read_segments.append((ref_start, ref_start + length - 1))
                ref_start += length
                query_start += length
            elif st == 2:
                # deletion in the read
                sj_ref.append((ref_start - 1, ref_start + length, 'deletion'))
                ref_start += length
            elif st == 1:
                # insertion
                sj_ref.append((ref_start - 1, ref_start - 1, 'insertion'))
                query_start += length
            elif st == 4:
                # soft-clipping
                query_start += length
        assert(len(read_segments) == len(sj_ref) + 1)
        final_read_segments = []
        final_sj_junctions = []

        for idx, jump in enumerate(sj_ref):
            left_segment = read_segments[idx]
            right_segment = read_segments[idx + 1]

            left_condition = (left_segment[1] - left_segment[0] >= 6)
            right_condition = (right_segment[1] - right_segment[0] >= 6)

            if left_condition:
                final_read_segments.append(left_segment)
            if right_condition:
                final_read_segments.append(right_segment)
            if left_condition and right_condition and jump[-1] == 'splice':
                final_sj_junctions.append(jump)

        return final_read_segments, final_sj_junctions, ref_start, read.query_length

    def assign_read_sgRNA(self, read) -> Tuple[str, Tuple[int, int]]:
        """Assign a read to a category of sgRNA.

        Args:
            read (pysam.AlignedSegment): an aligned read 

        Returns:
            str: a category: an ORF, non-canonical, or None (for non splicing reads).
            junction: splice junction for the assignment
        """
        if 'N' not in read.cigarstring:
            if read.get_overlap(*self.ORF1ab_start) > 0:
                return "1ab", None
            else:
                return None, None
        else:
            sjs, read_idx = self.get_spliced_junctions(read)
            for (_5prime, _3prime), (_, read_idx_3p) in zip(sjs, read_idx):
                query_first_ATG = read.query_sequence[read_idx_3p:].upper().find(
                    "ATG") + read_idx_3p
                if query_first_ATG == read_idx_3p-1:
                    return 'terminated', (_5prime, _3prime)
                ref_first_ATG = self.convert_query_index_to_ref(
                    query_first_ATG, read)
                if ref_first_ATG is None:
                    return 'unaligned', (_5prime, _3prime)
                if self.LeaderStart <= _5prime < self.LeaderEnd:
                    for sg_name, pos in self.sgRNAs.items():
                        if pos == ref_first_ATG:
                            return sg_name, (_5prime, _3prime)
            return 'non-canonical', (_5prime, _3prime)

    def count_sj(self, bam_file) -> Dict[str, int]:
        bam = pysam.AlignmentFile(bam_file, "rb")
        bin_sgRNA = {v: 0 for v in self.sgRNAs.keys()}
        bin_sgRNA['non-canonical'] = 0
        bin_sgRNA['1ab'] = 0
        bin_sgRNA['terminated'] = 0
        bin_sgRNA['unaligned'] = 0
        for read in bam.fetch(self.contig, self.LeaderStart, self.LeaderEnd):
            assignment, sj = self.assign_read_sgRNA(read)
            if assignment is not None:
                bin_sgRNA[assignment] += 1
        self.count_1ab(bam_file)
        return bin_sgRNA

    def count_1ab(self, bam_file) -> None:
        bam = pysam.AlignmentFile(bam_file, "rb")
        count_5prime = Counter()
        count_3prime = Counter()
        count = 0
        for read in bam.fetch(self.contig, 0, self.refLength):
            sjs, _ = self.get_spliced_junctions(read)
            if len(sjs) > 0:
                count_5prime.update([p5 for p5, _ in sjs])
                count_3prime.update([p3 for _, p3 in sjs])
            else:
                count += 1
        print("5prime counts --")
        for k, v in count_5prime.most_common(30):
            print(k, ': ', v)
        print("3prime counts --")
        for k, v in count_3prime.most_common(30):
            print(k, ': ', v)
        # for k in sorted(count_5prime.values()):
        #     print(f"{k}: {count_5prime[k]}")
        # print(f"ORF 1ab count is {count}")

    def count_sj_pos(self, pos, bam_file, max_depth=2147483647, paired=False) -> Counter:
        bam = pysam.AlignmentFile(bam_file, "rb")
        bin_sgRNA = []
        if paired:
            for r1, r2 in read_pair_generator(bam, f"{self.contig}:{pos-600}-{pos+600}", pos-1):
                allele = r1.query_sequence[self.convert_ref_index_to_query(
                    pos-1, r1)]
                if 'N' in r1.cigarstring:
                    assignment, sj = self.assign_read_sgRNA(r1)
                elif 'N' in r2.cigarstring:
                    assignment, sj = self.assign_read_sgRNA(r2)
                else:
                    assignment = None
                if assignment is not None:
                    bin_sgRNA.append((assignment, allele))
                else:
                    bin_sgRNA.append(("not_splice_junction", allele))
        else:
            for pileupColumn in bam.pileup(self.contig, pos-1, pos,
                                           stepper="samtools",
                                           min_base_quality=self.min_base_quality,
                                           min_mapping_quality=self.min_mapping_quality,
                                           max_depth=max_depth,
                                           ignore_overlaps=True,
                                           ignore_orhphans=True,
                                           truncate=True):
                for pileupRead in pileupColumn.pileups:
                    if not pileupRead.is_del and not pileupRead.alignment.is_duplicate:
                        if type(pileupRead.query_position) is int:
                            allele = pileupRead.alignment.query_sequence[pileupRead.query_position]
                        else:
                            allele = '_'
                        assignment, sj = self.assign_read_sgRNA(
                            pileupRead.alignment)
                        if assignment is not None:
                            bin_sgRNA.append((assignment, allele))
                        else:
                            bin_sgRNA.append(("not_splice_junction", allele))
        return Counter(bin_sgRNA)

    def get_allele_counts(self, bam_file, pos, max_depth=2147483647) -> Dict[str, int]:

        bam = pysam.AlignmentFile(bam_file, 'rb')
        allele_counts = {}

        for pileupColumn in bam.pileup(self.contig, pos-1, pos,
                                       stepper="samtools",
                                       min_base_quality=self.min_base_quality,
                                       min_mapping_quality=self.min_mapping_quality,
                                       max_depth=max_depth,
                                       ignore_overlaps=True,
                                       ignore_orhphans=True,
                                       truncate=True):
            for pileupRead in pileupColumn.pileups:
                if not pileupRead.is_del and not pileupRead.alignment.is_duplicate:
                    if type(pileupRead.query_position) is int:
                        allele = pileupRead.alignment.query_sequence[pileupRead.query_position]
                        if allele in allele_counts.keys():
                            allele_counts[allele] += 1
                        else:
                            allele_counts[allele] = 1

        return allele_counts

    def get_sj_reads(self, bam_file, mutList=[], paired=False) -> pd.DataFrame:
        bam = pysam.AlignmentFile(bam_file, "rb")
        count_5prime_3prime_pair = Counter()
        count = 0
#        if paired:
#            for r1, r2 in read_pair_generator(bam, f"{self.contig}:1-{self.refLength}"):
#                #sjs1, ref_pos1 = self.get_spliced_junctions_muts(r1, mutList)
#                #sjs2, _ = self.get_spliced_junctions_muts(r2, mutList)
#                _, sjs1, ref_pos1, _ = self.get_read_segments_junctions(r1)
#                _, sjs2, _, _ = self.get_read_segments_junctions(r2)
#                strand = not r1.is_reverse
#                if len(sjs1) > 0:
#                    count_5prime_3prime_pair.update([(p5, p3, strand) for p5, p3, entryType in sjs1
#                                                     if entryType == 'splice'])
#                if len(sjs2) > 0:
#                    count_5prime_3prime_pair.update([(p5, p3, strand) for p5, p3, entryType in sjs2
#                                                     if entryType == 'splice'])
#        else:
#            for read in bam.fetch(self.contig, 0, self.refLength):
#                if read.mapping_quality >= self.min_mapping_quality:
#                    sjs, _ = self.get_spliced_junctions_muts(read, mutList)
#                    #sjs, _ = self.get_spliced_junctions(read)
#                    strand = not read.is_reverse
#                    if len(sjs) > 0:
#                        count += 1
#                        count_5prime_3prime_pair.update([(p5, p3, strand)
#                                                         for p5, p3, entryType in sjs
#                                                         if entryType == 'splice'])

        print(f"{self.contig}\t{self.refLength}")
        for read in bam.fetch(self.contig, 0, self.refLength):
            if read.mapping_quality >= self.min_mapping_quality:
                #sjs, _ = self.get_spliced_junctions_muts(read, mutList)
                _, sjs, _, _ = self.get_read_segments_junctions(read)
                strand = not read.is_reverse
                if len(sjs) > 0:
                    count += 1
                    count_5prime_3prime_pair.update([(p5, p3, strand) for p5, p3, entryType in sjs if entryType == 'splice'])
        print(f"count is {count}")
        data = []
        for (p5, p3, strand), count in count_5prime_3prime_pair.items():
            data.append([p5, p3, strand, count])
        print(f"here")
        df = pd.DataFrame(data, columns=["5prime", "3prime", "strand", "count"])
        df["sgRNA_ref"] = df["3prime"].apply(lambda x: self.identify_ORF(self.reference.fetch(self.contig, x, self.refLength).find("ATG") + x))
        return df

    def get_phasing_reads(self, bam_file, splice_edges, width=0, paired=False):
        bam = pysam.AlignmentFile(bam_file, 'rb')
        phasing_reads = Counter()
        #phasing_reads = {}
        reject_count = 0
        print(f"ref legnth is {self.refLength}")
        if paired:
            #mean_read_len = IncrementalAverage()
            for reads in read_pair_generator(bam, f"{self.contig}:1-{self.refLength}"):
                s_plus = set()
                s_minus = set()

                ref_pos = None
                for read in reads:
                    read_seg, sjs, curr_ref_pos = self.get_read_segments(read)
                    #read_seg, sjs, curr_ref_pos, read_len = self.get_read_segments_junctions(read)
                    # mean_read_len.update(read_len)

                    if ref_pos and (read.reference_start > ref_pos):
                        break
                    else:
                        ref_pos = curr_ref_pos

                    for (left, right) in splice_edges:
                        for (_5prime, _3prime) in read_seg:
                            if (left < _3prime - width) and (right > _5prime):
                                s_minus.add((left, right))
                                break

                        for (_5prime, _3prime, _) in sjs:
                            if ((abs(left + float(width)/2 - _5prime) <= float(width)/2) and
                                    (abs(right + float(width)/2 - _3prime) <= float(width)/2)):
                                s_plus.add((left, right))
                                break

                if len(s_plus.intersection(s_minus)) == 0:
                    phasing_reads[tuple(
                        [tuple(sorted(s_plus)), tuple(sorted(s_minus))])] += 1
                else:
                    reject_count += 1
        else:
            for read in bam.fetch(self.contig, 0, self.refLength):
                if read.mapping_quality < self.min_mapping_quality:
                    reject_count += 1
                    continue

                s_plus = set()
                s_minus = set()

                read_seg, sjs, _ = self.get_read_segments(read)

                for (left, right) in splice_edges:
                    for (_5prime, _3prime) in read_seg:
                        if (left < _3prime - width) and (right > _5prime):
                            s_minus.add((left, right))
                            break

                    for (_5prime, _3prime, _) in sjs:
                        if ((abs(left + float(width)/2 - _5prime) <= float(width)/2) and (abs(right + float(width)/2 - _3prime) <= float(width)/2)):
                            s_plus.add((left, right))
                            break

                if len(s_plus.intersection(s_minus)) == 0:
                    phasing_reads[tuple([tuple(sorted(s_plus)), tuple(sorted(s_minus))])] += 1
                else:
                    reject_count += 1

        print(f"rejected {reject_count} reads")
        #print(f"average read length: {mean_read_len.mean}")
        # return dict(count_phasing_reads)
        return dict(phasing_reads)

    def get_sj_reads_phasing(self, bam_file, mutList=[], paired=False) -> Tuple[pd.DataFrame, List[List[Tuple[int, int, str]]]]:
        bam = pysam.AlignmentFile(bam_file, "rb")
        count_5prime_3prime_pair = Counter()
        phasing = []
        if paired:
            for r1, r2 in read_pair_generator(bam, f"{self.contig}:1-{self.refLength}"):
                sjs1, ref_pos1 = self.get_spliced_junctions_muts(r1, mutList)
                sjs2, _ = self.get_spliced_junctions_muts(r2, mutList)
                strand = not r1.is_reverse
                if len(sjs1) > 0:
                    count_5prime_3prime_pair.update([(p5, p3, strand) for p5, p3, entryType in sjs1
                                                     if entryType == 'splice'])
                if len(sjs2) > 0:
                    count_5prime_3prime_pair.update([(p5, p3, strand) for p5, p3, entryType in sjs2
                                                     if entryType == 'splice'])
                ref_pos2 = r2.reference_start
                if ref_pos2 > ref_pos1:
                    phasing.append(
                        sjs1 + [(ref_pos1, ref_pos2, 'pair')] + sjs2)
                else:
                    if len(sjs1) > 1:
                        phasing.append(sjs1)
                    if len(sjs2) > 1:
                        phasing.append(sjs2)
        else:
            for read in bam.fetch(self.contig, 0, self.refLength):
                if read.mapping_quality >= self.min_mapping_quality:
                    sjs, _ = self.get_spliced_junctions_muts(read, mutList)
                    #sjs, _ = self.get_spliced_junctions(read)
                    strand = not read.is_reverse
                    if len(sjs) > 0:
                        count_5prime_3prime_pair.update([(p5, p3, strand)
                                                         for p5, p3, entryType in sjs
                                                         if entryType == 'splice'])
                        #count_5prime_3prime_pair.update([(p5, p3, strand) for p5, p3 in sjs])
                    if len(sjs) > 1:
                        phasing.append(sjs)
        data = []
        for (p5, p3, strand), count in count_5prime_3prime_pair.items():
            data.append([p5, p3, strand, count])
        df = pd.DataFrame(
            data, columns=["5prime", "3prime", "strand", "count"])
        df["sgRNA_ref"] = df["3prime"].apply(lambda x: self.identify_ORF(
            self.reference.fetch(self.contig, x, self.refLength).find("ATG") + x))
        return df, phasing

    def save_all_sj_reads(self, bam_file, output_fname):
        df = self.get_sj_reads(bam_file)
        df.to_csv(output_fname, sep='\t', index=False)

    def save_all_sj_reads_phasing(self, bam_file, output_fname, output_phase):
        df, phasing = self.get_sj_reads_phasing(bam_file)
        df.to_csv(output_fname, sep='\t', index=False)

        with open(output_phase, 'w') as output:
            for idx, phase in enumerate(phasing):
                output.write(f">Phase_{idx}\n")
                for entry in phase:
                    output.write(f"{entry}\n")

    def save_sj_mut_spanning_count(self, bam_file, vcf_file, output_fname):
        vcf = pysam.VariantFile(vcf_file)
        variants = vcf.fetch()
        snv_pos = []
        allele_list = []
        orf_list = []
        count_list = []
        for var in variants:
            # is SNV
            if len(var.ref) == 1 and all(len(alt) == 1 for alt in var.alts) and var.qual >= self.min_snv_qual:
                bin_sgRNA = self.count_sj_pos(var.pos, bam_file)
                for (orf, allele), read_count in bin_sgRNA.items():
                    snv_pos.append(var.pos)
                    orf_list.append(orf)
                    allele_list.append(allele)
                    count_list.append(read_count)

        df = pd.DataFrame({'position': snv_pos,
                           'category': orf_list,
                           'allele': allele_list,
                           'count': count_list})

        df.to_csv(output_fname, sep='\t', index=False)


def main(args):
    SJcounter = SpliceJunction(
        args.fasta, min_base_qual=args.min_base_qual, min_mapping_qual=args.min_mapping_qual)
    if args.pos is None:
        bin_sgRNA = SJcounter.count_sj(args.bam)
    else:
        bin_sgRNA = SJcounter.count_sj_pos(args.pos, args.bam, args.pair_end)
    print(bin_sgRNA)


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bam", type=str, help="aligned bam file")
    parser.add_argument("-v", "--vcf", type=str, help="input vcf file")
    parser.add_argument("-f", "--fasta", type=str,
                        help="fasta file", required=True)
    parser.add_argument("--min-base-qual", type=int,
                        help="minimum base quality [20]", default=20)
    parser.add_argument("--min-mapping-qual", type=int,
                        help="minimum mapping quality [20]", default=20)
    parser.add_argument("--min-snv-qual", type=int,
                        help="minimum snv quality [20]", default=20)
    parser.add_argument("-p", "--pos", type=int,
                        help="position of the mutation, 1-based")
    parser.add_argument("-o", "--outputCSV", type=str, help="output csv file")
    parser.add_argument("--outputPhase", type=str,
                        help="output file for phasing reads")
    parser.add_argument("--paired", type=str2bool, nargs='?',
                        const=True, default=False, help="set to enable pair-end mode")
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    if not args.vcf is None and args.outputCSV is None:
        print("outputCSV is required with a vcf file input!")
        exit(1)
    elif not args.outputCSV is None and args.vcf is None:
        SJcounter = SpliceJunction(args.fasta, min_base_qual=args.min_base_qual,
                                   min_mapping_qual=args.min_mapping_qual, min_snv_qual=args.min_snv_qual)
        if args.outputPhase:
            SJcounter.save_all_sj_reads_phasing(
                args.bam, args.outputCSV, args.outputPhase)
        else:
            SJcounter.save_all_sj_reads(args.bam, args.outputCSV)
    elif args.outputCSV is None:
        main(args)
    else:
        SJcounter = SpliceJunction(args.fasta, min_base_qual=args.min_base_qual,
                                   min_mapping_qual=args.min_mapping_qual, min_snv_qual=args.min_snv_qual)
        SJcounter.save_sj_mut_spanning_count(
            args.bam, args.vcf, args.outputCSV)
