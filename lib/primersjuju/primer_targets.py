"""Primer targets specification from user.  Including Parsing and validated
primer_targets input file

"""
import re
from pycbio.tsv import TsvReader
from pycbio.hgdata.coords import Coords
from . import PrimersJuJuUserError

# for validation
REGION_COLS = frozenset(["region_5p", "region_3p"])
TRANSCRIPT_COLS = frozenset(["trans_track", "trans_id"])
REQUIRED_COLS = frozenset(frozenset(["target_id"]) | REGION_COLS | TRANSCRIPT_COLS)

class TargetTranscript:
    """a target transcript"""
    def __init__(self, trans_track, trans_id, user_attrs):
        self.trans_track = trans_track
        self.trans_id = trans_id
        self.user_attrs = user_attrs

class PrimerTarget:
    """
    A given primer target regions and associate target transcripts
    """
    def __init__(self, target_id, region_5p, region_3p, user_attrs):
        self.target_id = target_id
        self.region_5p = region_5p
        self.region_3p = region_3p
        self.user_attrs = user_attrs
        self.transcripts = {}  # by (track, id)

    def add_transcript(self, trans_track, trans_id, user_attrs):
        key = (trans_track, trans_id)
        if key in self.transcripts:
            raise PrimersJuJuUserError(f"duplicate transcript for primer: '{trans_track}', '{trans_id}'")
        trans = TargetTranscript(trans_track, trans_id, user_attrs)
        self.transcripts[key] = trans
        return trans

class PrimerTargets:
    """ All specified primer targets, each with a unique id"""
    def __init__(self):
        self.targets = {}

    def add_target(self, target_id, region_5p, region_3p, user_attrs):
        if target_id in self.targets:
            raise PrimersJuJuUserError(f"duplicate primer target_id '{target_id}")
        target = PrimerTarget(target_id, region_5p, region_3p, user_attrs)
        self.targets[target_id] = target
        return target

    def get_target(self, target_id):
        "target or None"
        return self.targets.get(target_id)

    def access_target(self, target_id):
        "target or error"
        target = self.targets.get(target_id)
        if target is None:
            raise PrimersJuJuUserError(f"unknown primer target_id references '{target_id}")
        return target

def _check_target_id(target_id):
    if not re.match("^[A-Za-z][-_.=%+A-Za-z0-9]*$", target_id):
        raise PrimersJuJuUserError(f"invalid target_id string '{target_id}', see documentation")

def _must_not_be_empty(columns, row):
    for col in columns:
        if row[col] == "":
            raise PrimersJuJuUserError(f"row column {col} must not be empty")

def _must_be_empty(columns, row):
    for col in columns:
        if row[col] != "":
            raise PrimersJuJuUserError(f"row column {col} must be empty")

def _parse_coords(coord_str):
    coords = Coords.parse(coord_str, oneBased=True)
    if len(coords) > 10000000:
        raise PrimersJuJuUserError(f"coordinates seem absurdly long: '{coord_str}'")
    return coords

def _get_user_cols(rows):
    row0 = rows[0]
    target_user_cols = set()
    transcript_user_cols = set()
    for col in row0._columns_:
        if col not in REQUIRED_COLS:
            if col.startswith("trans_"):
                transcript_user_cols.add(col)
            else:
                target_user_cols.add(col)
    return target_user_cols, transcript_user_cols

def _build_column_dict(columns, row):
    return {col: row[col] for col in columns}

def _do_add_primary_row(primer_targets, target_user_cols, transcript_user_cols, row):
    _check_target_id(row.target_id)
    _must_not_be_empty(REQUIRED_COLS, row)

    target = primer_targets.add_target(row.target_id, row.region_5p, row.region_3p,
                                       _build_column_dict(target_user_cols, row))
    target.add_transcript(row.trans_track, row.trans_id,
                          _build_column_dict(transcript_user_cols, row))

def _add_primary_row(primer_targets, target_user_cols, transcript_user_cols, row):
    try:
        _do_add_primary_row(primer_targets, target_user_cols, transcript_user_cols, row)
    except Exception as ex:
        raise PrimersJuJuUserError(f"error parsing primary row: {str(row)}") from ex

def _do_add_continue_row(primer_targets, transcript_user_cols, row):
    _check_target_id(row.target_id)
    _must_be_empty(REGION_COLS, row)
    _must_not_be_empty(TRANSCRIPT_COLS, row)
    target = primer_targets.access_target(row.target_id)
    target.add_transcript(row.trans_track, row.trans_id,
                          _build_column_dict(transcript_user_cols, row))

def _add_continue_row(primer_targets, transcript_user_cols, row):
    try:
        _do_add_continue_row(primer_targets, transcript_user_cols, row)
    except Exception as ex:
        raise PrimersJuJuUserError(f"error parsing continuation row: {str(row)}") from ex

def _primer_targets_build(rows):
    # since no order required, two passes; one to for primary data and the
    # other for continuation
    primer_targets = PrimerTargets()
    target_user_cols, transcript_user_cols = _get_user_cols(rows)
    for row in rows:
        if row.region_5p != "":
            _add_primary_row(primer_targets, target_user_cols, transcript_user_cols, row)
    for row in rows:
        if row.region_5p == "":
            _add_continue_row(primer_targets, transcript_user_cols, row)

def _check_required_columns(rows):
    if len(rows) == 0:
        raise PrimersJuJuUserError("no data in TSV")
    row0 = rows[0]
    for col in REQUIRED_COLS:
        if col not in row0:
            raise PrimersJuJuUserError(f"required column is missing: '{col}'")

def primer_targets_read(primer_targets_tsv):
    """read all primer targets into PrimerTargets object"""
    try:
        rows = [row for row in TsvReader(primer_targets_tsv)]
        _check_required_columns(rows)
        return _primer_targets_build(rows)
    except Exception as ex:
        raise PrimersJuJuUserError(f"error parsing primary target specification TSV: '{primer_targets_tsv}'") from ex