import pandas as pd
import pysam
from bam_io import SpliceJunction as SJ
import argparse
import sys
from typing import List, Dict, Tuple

"""
!!Alert: Always use [start, end) intervals, and 0-based coordinate internally
"""


class Exon:
    """
    A exon in a transcript
    """

    def __init__(self, index, start, end):
        self.index = index
        self.start = start
        self.end = end


class Transcript:
    """
    A transcript
    """

    def __init__(self, name, start, end):
        self.name = name
        self.start = start
        self.end = end
        self.exons = []

    def add_exon(self, index, start, end):
        self.exons.append(Exon(index, start, end))

    def assign_read(self, read):
        pass

class Segments:
    def __init__(self, transcriptome: Dict[str, Transcript], length):
        set_breakpoints = set()
        for transcript in transcriptome.values():
            for exon in transcript.exons:
                set_breakpoints.update([exon.start, exon.end])
        self.breakpoints = sorted(list(set_breakpoints))
        self.n_segments = 1 + len(self.breakpoints)
        
        if self.breakpoints[0] > 0:
            self.segments = [(0, self.breakpoints[0])]
        else:
            self.segments = []
        self.segments += [(v, w) for v, w in zip(self.breakpoints[:-1], self.breakpoints[1:])]
        if self.breakpoints[-1] < length:
            self.segments.append((self.breakpoints[-1], length))
        print('\n'.join(list(f"{x}: {x[1] - x[0]}" for x in self.segments)))


def read_GTF(transcript_file) -> pd.DataFrame:
    with open(transcript_file) as ifile:
        data_transcript = []
        header_names = ['seqname', 'source', 'feature',
                        'start', 'end', 'score', 'strand', 'frame']
        for line in ifile:
            record = {k: v for k, v in zip(header_names, line.split('\t')[:8])}
            # Use 0-based coordinate
            if record["start"] is not None and len(record["start"]) > 0:
                record["start"] = int(record["start"]) - 1
            record["end"] = int(record["end"])
            for x in line.split('\t')[8].split(";"):
                if len(x.strip()) > 0:
                    record.update({x.split()[0].strip(): x.split()[1].strip()})
            data_transcript.append(record)
        return pd.DataFrame(data=data_transcript)


def build_transcript(df_transcript: pd.DataFrame) -> Dict[str, Transcript]:
    dict_transcript = {}
    for _, row in df_transcript.loc[df_transcript["feature"] == "transcript"].iterrows():
        dict_transcript[row["transcript_id"]] = Transcript(row["transcript_id"], row["start"], row["end"])
    for _, row in df_transcript.loc[df_transcript["feature"] == "exon"].iterrows():
        dict_transcript[row["transcript_id"]].add_exon(row["exon"], row["start"], row["end"])
    return dict_transcript


def read_bam(bam_file, contig_name, ref_file):
    SJcounter = bam_io.SpliceJunction(ref_file, contig = contig_name)
    df_sj_reads = SJcounter.get_sj_reads(bam_file, paired=False)


def get_read_segments_junctions(read):
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
                sj_ref.append((ref_start, ref_start + length - 1, "segment"))
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
        final_read_segments = []
        final_sj_junctions = []
        
        for idx, jump in enumerate(sj_ref):
            if jump[-1] == "splice":
                left_segment = sj_ref[idx-1]
                right_segment = sj_ref[idx + 1]
                
                left_condition = (left_segment[1] - left_segment[0] >= 6)
                right_condition = (right_segment[1] - right_segment[0] >= 6)
            
                if left_condition:
                    final_read_segments.append(left_segment)
                if right_condition:
                    final_read_segments.append(right_segment)
                if left_condition and right_condition and jump[-1] == 'splice':
                    final_sj_junctions.append(jump)
        
        return final_read_segments, final_sj_junctions, ref_start, read.query_length


def main(transcript_file, bam_file, contig="chrSCV"):
    # df_transcript = read_GTF(transcript_file)
    # transcriptome = build_transcript(df_transcript)
    # segments = Segments(transcriptome, 29903)

    transcriptome = {}
    with open(transcript_file) as ifile:
        name = None
        for line in ifile:
            if line.startswith(">"):
                if name is None:
                    data_transcript = []
                else:
                    transcriptome[name] = pd.DataFrame(data_transcript)
                    data_transcript = []
                name = line[1:]
            else:
                former, latter, type_ = line.split()
                former = former.strip('][').split(',')
                latter = latter.strip('][').split(',')
                data_transcript.append({
                    "v_start": former[0],
                    "v_end": former[1],
                    "w_start": latter[0],
                    "w_end": latter[1],
                    "type": type_})

    linked_splices = {}
    for name, transcript in transcriptome.items():
        if (transcript["type"] == "splice").sum() > 1:
            splices = transcript.loc[transcript["type"] == "splice", ["v_end", "w_start"]].astype(int).to_numpy()
            linked_splices[name] = splices
            print(name, transcript[transcript["type"] == "splice"], splices)

    tol = 10

    bam = pysam.AlignmentFile(bam_file)
    counts_transcript = {name: {"exact": 0, "partial": 0} for name in linked_splices.keys()}
    cnt = 0
    for read in bam.fetch(contig):
        cnt += 1
        if cnt % 1000 == 0:
            print(f"{cnt}", end='\r')
        if read.cigarstring.count("N") >= 2 and read.mapping_quality > 30:
            _, jumps, _, _ = get_read_segments_junctions(read)
            for transcript, splices in linked_splices.items():
                for start, end, kind in jumps:
                    if kind == "splice":
                        for row in splices:
                            if abs(start - row[0]) < tol and abs(end - row[1]) < tol:
                                is_transcript_failed = False
                                break
                        else:
                            is_transcript_failed = True
                    if is_transcript_failed:
                        # hits no splice
                        break
                else:
                    # all jumps succeed
                    if read.cigarstring.count("N") == splices.shape[0]:
                        counts_transcript[transcript]["exact"] += 1
                    else:
                        counts_transcript[transcript]["partial"] += 1
    print(counts_transcript)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Verify inferred transcripts with long reads.")
    parser.add_argument("-g", "--transcript", type=str,
                        help="input GTF containing transcripts")
    parser.add_argument("-b", "--bam", type=str, help="BAM file of long reads")
    parser.add_argument("-r", "--ref", type=str, help="FASTA file of reference genome, assume only containing SARS-CoV-2")
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args.transcript, args.bam)
