# Class for configuration objects in the python code based configuration file
# your config.py should import this and construct a instance of
# PrimersJuJuConfig in a variable named 'config'
#
from dataclasses import dataclass
from typing import Optional, Tuple
import primer3
from primersjuju.uniqueness_query import UniquenessQuery, IsPcrServerSpec
from primersjuju.genome_data import GenomeData

class GenomeConfig:
    """Configuration for a particular assembly.  Setting either of ispcr specs to None
    cases the corresponding off target query to be skipped."""

    def __init__(self,
                 genome_data: GenomeData,
                 genome_ispcr_spec: IsPcrServerSpec,
                 transcriptome_ispcr_spec: IsPcrServerSpec):
        self.genome_data = genome_data
        self.genome_ispcr_spec = genome_ispcr_spec
        self.transcriptome_ispcr_spec = transcriptome_ispcr_spec
        self.__uniqueness_query = None   # lazy

    @property
    def genome_name(self):
        return self.genome_data.genome_name

    @property
    def uniqueness_query(self) -> UniquenessQuery:
        "lazy get/create UniquenessQuery object, or None if not configured"
        if ((self.__uniqueness_query is None) and
            ((self.genome_ispcr_spec is not None) or (self.transcriptome_ispcr_spec is not None))):
            self.__uniqueness_query = UniquenessQuery(self.genome_data,
                                                      self.genome_ispcr_spec, self.transcriptome_ispcr_spec)
        return self.__uniqueness_query

class Primer3ThermoArgs(primer3.argdefaults.Primer3PyArguments):
    """Primer3 thermoanalysis arguments with standard defaults."""
    pass

@dataclass
class Primer3ThermoFitlers:
    """Primer3 thermoanalysis filters.
    For floating range tests, values are either tuples of (min, max) or None to not run the test.
    For stability tests, it is the minimum delta-G or None
    """
    tm_range: Optional[Tuple[float, float]] = None
    hairpin_min_dg: Optional[float] = None
    homodimer_min_dg: Optional[float] = None
    heterodimer_min_dg: Optional[float] = None
    end_stability_min_dg: Optional[float] = None

class DefaultThermoFitlers(Primer3ThermoFitlers):
    """A default set of Primer3ThermoFitlers"""
    def __init__(self):
        super().__init__(tm_range=(55.0, 70.0),
                         hairpin_min_dg=-3.0,
                         homodimer_min_dg=-3.0,
                         heterodimer_min_dg=-5.0,
                         end_stability_min_dg=-10.0)


class Primer3Config:
    """Options to pass to primer3, see primer3 manual Global Input Tags for a
    description.
    https://primer3.org/manual.html#globalTags
    """
    def __init__(self):
        self.num_5_prime_strong_match = 0   # builds PRIMER_MUST_MATCH_FIVE_PRIME
        self.PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION = 8
        self.PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION = 8
        self.PRIMER_MAX_POLY_X = None
        self.PRIMER_MAX_SELF_ANY = None
        self.PRIMER_MAX_SELF_ANY_TH = None
        self.PRIMER_PAIR_MAX_COMPL_ANY = None
        self.PRIMER_PAIR_MAX_COMPL_ANY_TH = None
        self.PRIMER_PICK_INTERNAL_OLIGO = 0
        self.PRIMER_OPT_SIZE = 20
        self.PRIMER_MIN_SIZE = 18
        self.PRIMER_MAX_SIZE = 22
        self.PRIMER_NUM_RETURN = 5

        # library files
        self.misprime_lib = None
        self.mishyb_lib = None

        ##
        # thermoanalysis
        ##
        self.thermo_filters = None  # Primer3ThermoFilters, set of filter by thermodynamic attributes
        self.thermo_args = None  # Primer3ThermoArgs, use defaults if None


_default_primer3_config = Primer3Config()

class PrimersJuJuConfig:
    """Main configuration object. An instance of this object must be create
    and stored in a variable Config"""
    def __init__(self, *, primer3_config=_default_primer3_config):
        self.genomes = {}
        self.primer3 = primer3_config

        # derived, not set in config file
        self.genome = None  # selected genome, set at load from genomes

    def add_genome(self, genome_config: GenomeConfig):
        self.genomes[genome_config.genome_name] = genome_config
