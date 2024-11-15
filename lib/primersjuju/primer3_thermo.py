"""
Interface to primer3 thermodynamics functions exported by primer3_py
"""
from collections import namedtuple
from dataclasses import dataclass
from typing import Optional
import primer3
from .config import Primer3ThermoArgs


class ThermoEval(namedtuple("ThermoEval", ("passed", "val"))):
    """Result for thermodynamic checks.  The passed field is True if the check
    was within the configured range.  The val field are either a TM or a
    delta-G.
    """
    __slots__ = ()

class ThermoEvalPair(namedtuple("ThermoEvalPair", ("fwd", "rev"))):
    """Forward and reverse values thermodynamic check results"""

    @property
    def passed(self):
        "did both pass?"
        return self.fwd.passed and self.rev.passed

def passed(teval):
    return (teval is None) or teval.passed

@dataclass
class Primer3ThermoResults:
    """
    Results from primer3 thermodynamics checks.  Fields are None if the check
    was not done, or contain a pass field, plus a value.
    """

    tm: Optional[ThermoEvalPair] = None
    hairpin: Optional[ThermoEvalPair] = None
    homodimer: Optional[ThermoEvalPair] = None
    heterodimer: Optional[ThermoEval] = None
    end_stability: Optional[ThermoEval] = None

    @property
    def passed(self):
        "do all checks pass or are None?"
        return (passed(self.tm) and passed(self.hairpin) and passed(self.homodimer) and
                passed(self.heterodimer) and passed(self.end_stability))


def eval_tm1(seq, args, filters):
    tm = primer3.calc_tm(seq,
                         mv_conc=args.mv_conc,
                         dv_conc=args.dv_conc,
                         dntp_conc=args.dntp_conc,
                         dna_conc=args.dna_conc,
                         dmso_conc=args.dmso_conc,
                         dmso_fact=args.dmso_fact,
                         formamide_conc=args.formamide_conc,
                         annealing_temp_c=args.annealing_temp_c,
                         max_nn_length=args.max_nn_length,
                         tm_method=args.tm_method,
                         salt_corrections_method=args.salt_corrections_method)
    return ThermoEval(filters.tm_range[0] <= tm <= filters.tm_range[1], tm)

def eval_tm(primer3_pair, args, filters):
    return ThermoEvalPair(eval_tm1(primer3_pair.PRIMER_LEFT_SEQUENCE, args, filters),
                          eval_tm1(primer3_pair.PRIMER_RIGHT_SEQUENCE, args, filters))

def eval_hairpin1(seq, args, filters):
    tr = primer3.calc_hairpin(seq,
                              mv_conc=args.mv_conc,
                              dv_conc=args.dv_conc,
                              dntp_conc=args.dntp_conc,
                              dna_conc=args.dna_conc,
                              temp_c=args.temp_c,
                              max_loop=args.max_loop,
                              output_structure=False)
    tr.check_exc()
    return ThermoEval(tr.dg >= filters.hairpin_min_dg, tr.dg)

def eval_hairpin(primer3_pair, args, filters):
    return ThermoEvalPair(eval_hairpin1(primer3_pair.PRIMER_LEFT_SEQUENCE, args, filters),
                          eval_hairpin1(primer3_pair.PRIMER_RIGHT_SEQUENCE, args, filters))

def eval_homodimer1(seq, args, filters):
    tr = primer3.calc_homodimer(seq,
                                mv_conc=args.mv_conc,
                                dv_conc=args.dv_conc,
                                dntp_conc=args.dntp_conc,
                                dna_conc=args.dna_conc,
                                temp_c=args.temp_c,
                                max_loop=args.max_loop,
                                output_structure=False)
    tr.check_exc()
    return ThermoEval(tr.dg >= filters.homodimer_min_dg, tr.dg)

def eval_homodimer(primer3_pair, args, filters):
    return ThermoEvalPair(eval_homodimer1(primer3_pair.PRIMER_LEFT_SEQUENCE, args, filters),
                          eval_homodimer1(primer3_pair.PRIMER_RIGHT_SEQUENCE, args, filters))

def eval_heterodimer(primer3_pair, args, filters):
    tr = primer3.calc_heterodimer(primer3_pair.PRIMER_LEFT_SEQUENCE,
                                  primer3_pair.PRIMER_RIGHT_SEQUENCE,
                                  mv_conc=args.mv_conc,
                                  dv_conc=args.dv_conc,
                                  dntp_conc=args.dntp_conc,
                                  dna_conc=args.dna_conc,
                                  temp_c=args.temp_c,
                                  max_loop=args.max_loop,
                                  output_structure=False)
    tr.check_exc()
    return ThermoEval(tr.dg >= filters.heterodimer_min_dg, tr.dg)

def eval_end_stability(primer3_pair, args, filters):
    tr = primer3.calc_end_stability(primer3_pair.PRIMER_LEFT_SEQUENCE,
                                    primer3_pair.PRIMER_RIGHT_SEQUENCE,
                                    mv_conc=args.mv_conc,
                                    dv_conc=args.dv_conc,
                                    dntp_conc=args.dntp_conc,
                                    dna_conc=args.dna_conc,
                                    temp_c=args.temp_c,
                                    max_loop=args.max_loop)
    tr.check_exc()
    return ThermoEval(tr.dg >= filters.end_stability_min_dg, tr.dg)

def eval_thermodynamics(primer3_pair, args, filters):
    results = Primer3ThermoResults()
    if filters.tm_range is not None:
        results.tm = eval_tm(primer3_pair, args, filters)
    if filters.hairpin_min_dg is not None:
        results.hairpin = eval_hairpin(primer3_pair, args, filters)
    if filters.homodimer_min_dg is not None:
        results.homodimer = eval_homodimer(primer3_pair, args, filters)
    if filters.heterodimer_min_dg is not None:
        results.heterodimer = eval_heterodimer(primer3_pair, args, filters)
    if filters.end_stability_min_dg is not None:
        results.end_stability = eval_end_stability(primer3_pair, args, filters)
    return results

def evalulate_thermodynamics(primer3_pair, args, filters):
    "evaluated thermodynamics and flag rejected if they fail"
    primer3_pair.thermo_results = eval_thermodynamics(primer3_pair, args, filters)
    primer3_pair.passed = primer3_pair.thermo_results.passed

def primer3_thermo_evals(primer3_results, primer3_config):
    args = primer3_config.thermo_args
    if args is None:
        args = Primer3ThermoArgs()
    for primer3_pair in primer3_results.pairs:
        evalulate_thermodynamics(primer3_pair, args, primer3_config.thermo_filters)
