# common output functions

import os.path as osp
# import urllib.parse
from pycbio.sys import fileOps
from pycbio.sys.svgcolors import SvgColors
from pycbio.hgdata.bed import Bed
from primersjuju.transcript_features import features_to_genomic_coords
from primersjuju.target_transcripts import ExonFeature
from primersjuju.primer3_interface import primer3_dump_args, primer3_annotate_rna

TARGETS_COLOR = SvgColors.blue
PRIMERS_ON_COLOR = SvgColors.darkmagenta
PRIMERS_OFF_GENOME_COLOR = SvgColors.orange
PRIMERS_OFF_TRANSCRIPTOME_COLOR = SvgColors.red
UNIQ_ON_COLOR = SvgColors.blue
UNIQ_OFF_COLOR = SvgColors.red
UNIQ_NON_COLOR = SvgColors.yellow


def _coords_to_bed(name, color, coords_list, extra_cols=None):
    coords_list = sorted(coords_list, key=lambda c: (c.name, c.start, c.end))
    first = coords_list[0]
    last = coords_list[-1]

    bed = Bed(first.name, first.start, last.end, name,
              strand=first.strand, thickStart=first.start, thickEnd=last.end,
              itemRgb=color.toRgb8Str(), numStdCols=12,
              extraCols=extra_cols)
    for coords in coords_list:
        bed.addBlock(coords.start, coords.end)
    return bed

def build_target_bed(target_transcripts):
    coords_list = [target_transcripts.region_5p, target_transcripts.region_3p]
    return _coords_to_bed(target_transcripts.target_id, TARGETS_COLOR, coords_list)

_primer_bed_columns = (
    'PRIMER_PAIR_PRODUCT_SIZE',    # 1133
    'PRIMER_PAIR_COMPL_ANY_TH',    # 12.498152433966595
    'PRIMER_PAIR_COMPL_END_TH',    # 5.033689592343535
    'PRIMER_PAIR_PENALTY',         # 0.42111006752037383
    'PRIMER_LEFT',                 # (27, 20)
    'PRIMER_LEFT_SEQUENCE',        # 'GGTTCTTCTGCGCTACTGCT'
    'PRIMER_LEFT_END_STABILITY',   # 4.24
    'PRIMER_LEFT_GC_PERCENT',      # 55.0
    'PRIMER_LEFT_HAIRPIN_TH',      # 36.97289132676252
    'PRIMER_LEFT_PENALTY',         # 0.3904716385882807
    'PRIMER_LEFT_SELF_ANY_TH',     # 9.732052967213463
    'PRIMER_LEFT_SELF_END_TH',     # 0.0
    'PRIMER_LEFT_TM',              # 60.39047163858828
    'PRIMER_RIGHT',                # (1159, 20)
    'PRIMER_RIGHT_SEQUENCE',       # 'CAAAAACCCACGCAGACAGG'
    'PRIMER_RIGHT_END_STABILITY',  # 4.0
    'PRIMER_RIGHT_GC_PERCENT',     # 55.0
    'PRIMER_RIGHT_HAIRPIN_TH',     # 0.0
    'PRIMER_RIGHT_PENALTY',        # 0.03063842893209312
    'PRIMER_RIGHT_SELF_ANY_TH',    # 0.0
    'PRIMER_RIGHT_SELF_END_TH',    # 0.0
    'PRIMER_RIGHT_TM',             # 59.96936157106791
)

def _get_extra_cols(primer_design):
    """get extra BED columns from primer design"""
    extra_cols = [primer_design.transcriptome_on_target_cnt(),
                  primer_design.transcriptome_off_target_cnt(),
                  primer_design.transcriptome_non_target_cnt(),
                  primer_design.genome_on_target_cnt(),
                  primer_design.genome_off_target_cnt(),
                  primer_design.genome_non_target_cnt()]
    # directly from primer3
    for col_name in _primer_bed_columns:
        col = getattr(primer_design.primer3_pair, col_name)
        if isinstance(col, float):
            col = "{:.2f}".format(col)
        extra_cols.append(col)
    return extra_cols

def _primer_color(primer_design):
    if primer_design.transcriptome_off_target_cnt() > 0:
        return PRIMERS_OFF_TRANSCRIPTOME_COLOR
    elif primer_design.genome_off_target_cnt() > 0:
        return PRIMERS_OFF_GENOME_COLOR
    else:
        return PRIMERS_ON_COLOR

def _primer_to_bed(primer_design):
    coords_list = (features_to_genomic_coords(primer_design.features_5p, ExonFeature) +
                   features_to_genomic_coords(primer_design.features_3p, ExonFeature))
    return _coords_to_bed(primer_design.ppair_id, _primer_color(primer_design), coords_list,
                          _get_extra_cols(primer_design))

def build_primer_beds(primer_designs):
    return sorted([_primer_to_bed(pd) for pd in primer_designs.designs],
                  key=Bed.genome_sort_key)

def _genome_hit_to_bed(hit, name, color):
    return _coords_to_bed(name, color, (hit.left_coords, hit.right_coords))

def _genome_hits_to_bed(hits, name, color):
    if hits is None:
        return []  # no data
    return [_genome_hit_to_bed(hit, name, color) for hit in hits]

def _transcriptome_hit_to_bed(hit, name_pre, color):
    return _coords_to_bed(name_pre + '|' + hit.trans_id + '|' + hit.gene_name, color,
                          features_to_genomic_coords(hit.left_features, ExonFeature) +
                          features_to_genomic_coords(hit.right_features, ExonFeature))

def _transcriptome_hits_to_bed(hits, name_pre, color):
    if hits is None:
        return []  # no data
    return [_transcriptome_hit_to_bed(hit, name_pre, color) for hit in hits]

def _build_design_uniqueness_hits_beds(primer_design):
    beds = []
    beds += _genome_hits_to_bed(primer_design.genome_on_targets, "on=" + primer_design.ppair_id, UNIQ_ON_COLOR)
    beds += _genome_hits_to_bed(primer_design.genome_off_targets, "off=" + primer_design.ppair_id, UNIQ_OFF_COLOR)
    beds += _genome_hits_to_bed(primer_design.genome_off_targets, "non=" + primer_design.ppair_id, UNIQ_NON_COLOR)
    beds += _transcriptome_hits_to_bed(primer_design.transcriptome_on_targets, "on=" + primer_design.ppair_id, UNIQ_ON_COLOR)
    beds += _transcriptome_hits_to_bed(primer_design.transcriptome_off_targets, "off=" + primer_design.ppair_id, UNIQ_OFF_COLOR)
    beds += _transcriptome_hits_to_bed(primer_design.transcriptome_non_targets, "non=" + primer_design.ppair_id, UNIQ_NON_COLOR)
    return beds

def build_uniqueness_hits_beds(primer_designs):
    beds = []
    for primer_design in primer_designs.designs:
        beds += _build_design_uniqueness_hits_beds(primer_design)
    return sorted(beds, key=Bed.genome_sort_key)

def _write_beds(beds, bed_file):
    with open(bed_file, "w") as fh:
        for bed in beds:
            bed.write(fh)

def _get_out_path(outdir, target_transcripts, suffix):
    return osp.join(outdir, target_transcripts.target_id + "." + suffix)

def output_target_design_file(outdir, target_transcripts):
    """get path to target design TSV file for an experiment.  The existence of this file indicates design
    for this target was complete.
    """
    return _get_out_path(outdir, target_transcripts, ".design.tsv")

_design_tsv_header = ("target_id", "design_status", "transcript_id",
                      "primer_id", "left_primer", "right_primer",
                      "transcriptome_on_targets", "transcriptome_off_targets",
                      "genome_on_targets", "genome_off_targets", "position",
                      "browser", "annotated_rna")

def _get_position(primer_designs):
    return str(primer_designs.target_transcript.trans_coords)

def _make_browser_link(primer_designs, hub_urls):
    browser_url = "https://genome.ucsc.edu/cgi-bin/hgTracks"

    cgi_args = ["position=" + _get_position(primer_designs),
                "db=" + primer_designs.target_transcripts.genome_name,
                "genome=" + primer_designs.target_transcripts.genome_name]
    for hub_url in hub_urls:
        cgi_args.append("hubUrl=" + hub_url)
    # FIXME: browser can't handle decoding
    # return browser_url + "?" + "&".join([urllib.parse.quote(a) for a in cgi_args])
    return browser_url + "?" + "&".join(cgi_args)

def _make_browser_excel_link(primer_designs, hub_urls):
    url = _make_browser_link(primer_designs, hub_urls)
    return f'=HYPERLINK("{url}", "{_get_position(primer_designs)}")'

def _write_primer_pair_design(fh, primer_designs, primer_design, first, hub_urls):
    "write one design, if primer_design is None, it means there were no primers found"
    row = [primer_designs.target_transcripts.target_id,
           primer_designs.status, primer_designs.target_transcripts.transcripts[0].trans_id]
    if primer_design is None:
        row += 7 * ['']
    else:
        row += [primer_design.ppair_id,
                primer_design.primer3_pair.PRIMER_LEFT_SEQUENCE,
                primer_design.primer3_pair.PRIMER_RIGHT_SEQUENCE,
                primer_design.transcriptome_on_target_cnt(),
                primer_design.transcriptome_off_target_cnt(),
                primer_design.genome_on_target_cnt(),
                primer_design.genome_off_target_cnt()]
    if first:
        row.append(_get_position(primer_designs))
        row.append(_make_browser_excel_link(primer_designs, hub_urls) if hub_urls is not None else "")
        row.append(primer3_annotate_rna(primer_designs.target_transcripts.transcripts[0]))
    else:
        row += 3 * ['']
    fileOps.prRow(fh, row)

def _write_primer_designs(primer_designs, tsv_file, hub_urls):
    with open(tsv_file, "w") as fh:
        fileOps.prRow(fh, _design_tsv_header)
        first = True
        if len(primer_designs.designs) == 0:
            _write_primer_pair_design(fh, primer_designs, None, first, hub_urls)
        else:
            for primer_design in primer_designs.designs:
                _write_primer_pair_design(fh, primer_designs, primer_design, first, hub_urls)
                first = False

def output_target_designs(outdir, target_transcripts, primer_designs, hub_urls=None):
    """output primer TSV, BEDs, and debug information for one target.  The
    $target_id.design.tsv file is created atomically, so it can be used as a
    marker that this design is complete.  If hub_urls is specified, the
    design.tsv file will contain links to the target region in the UCSC
    browser.
    """

    target_transcript = target_transcripts.transcripts[0]
    fileOps.ensureDir(outdir)

    # debugging
    with open(_get_out_path(outdir, target_transcripts, "debug.txt"), 'w') as fh:
        print("target_id", target_transcripts.target_id, file=fh)
        print("annotated_rna:", primer3_annotate_rna(target_transcripts.transcripts[0]), file=fh)
        primer3_dump_args(fh, target_transcript)
        print(file=fh)
        primer_designs.primer3_results.dump(fh)
        print(file=fh)
        primer_designs.dump(fh)
        print(file=fh)

    # BEDs
    _write_beds([build_target_bed(target_transcripts)],
                _get_out_path(outdir, target_transcripts, "target.bed"))
    _write_beds(build_primer_beds(primer_designs),
                _get_out_path(outdir, target_transcripts, "primers.bed"))
    _write_beds(build_uniqueness_hits_beds(primer_designs),
                _get_out_path(outdir, target_transcripts, "uniqueness.bed"))

    with fileOps.AtomicFileCreate(_get_out_path(outdir, target_transcripts, "designs.tsv")) as tmp_tsv:
        _write_primer_designs(primer_designs, tmp_tsv, hub_urls)