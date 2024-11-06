import pytest
import os.path as osp
import copy
from . import mydir
from primersjuju.config import PrimersJuJuConfig, DefaultThermoFitlers, Primer3ThermoArgs
from primersjuju.genome_data import GenomeData
from primersjuju.primer_target_spec import primer_target_specs_read
from primersjuju.uniqueness_query import IsPcrServerSpec, UniquenessQuery

def _test_data_file(fname):
    return osp.join(mydir, "../data", fname)

@pytest.fixture(scope="session")
def config_hg38():
    test_gdata = GenomeData("hg38",
                            _test_data_file("hg38.2bit"),
                            assembly_report=_test_data_file("GCF_000001405.40_GRCh38.p14_assembly_report.txt"))
    test_gdata.add_track("gencodeV39",
                         _test_data_file("gencodeV39.bb"),
                         "https://hgdownload.soe.ucsc.edu/gbdb/hg38/gencode/gencodeV39.bb")
    test_gdata.add_track("WTC11_consolidated",
                         _test_data_file("WTC11_consolidated.bigBed"),
                         "http://conesalab.org/LRGASP/LRGASP_hub/hg38/Human_samples/WTC11_consolidated.bigBed")
    config = PrimersJuJuConfig()
    config.add_genome(test_gdata)
    config.genome = test_gdata
    return config

@pytest.fixture(scope="session")
def config_hg38_thermo(config_hg38):
    config = copy.deepcopy(config_hg38)
    config.primer3.thermo_filters = DefaultThermoFitlers()
    config.primer3.thermo_args = Primer3ThermoArgs()
    return config

@pytest.fixture(scope="session")
def example_targets_specs():
    return primer_target_specs_read(osp.join(mydir, "../../docs/primer-targets-example.tsv"))

@pytest.fixture(scope="session")
def wtc11_track(config_hg38):
    return config_hg38.genome.get_track("WTC11_consolidated")

@pytest.fixture(scope="session")
def wtc11_targets_specs_set1():
    return primer_target_specs_read(osp.join(mydir, "../input/LRGASP_WTC11_Target_set1.tsv"))

@pytest.fixture(scope="session")
def hg38_ispcr_spec():
    # hg38
    return IsPcrServerSpec("blat1d.soe.ucsc.edu", 17903, osp.join(mydir, "../data"))

@pytest.fixture(scope="session")
def gencode_ispcr_spec():
    # hg38KgSeqV39
    return IsPcrServerSpec("blat1a.soe.ucsc.edu", 17907, osp.join(mydir, "../data"),
                           trans_bigbed=osp.join(mydir, "../data/gencodeV39.bb"))

@pytest.fixture(scope="session")
def hg38_uniqueness_query(config_hg38, hg38_ispcr_spec, gencode_ispcr_spec):
    return UniquenessQuery(config_hg38.genome, hg38_ispcr_spec, gencode_ispcr_spec)
