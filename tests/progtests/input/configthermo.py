# configuration file with hg38

from os import path as osp
from primersjuju.config import PrimersJuJuConfig, GenomeConfig, Primer3Config, Primer3ThermoArgs, DefaultThermoFitlers
from primersjuju.genome_data import GenomeData
from primersjuju.uniqueness_query import IsPcrServerSpec

data_dir = osp.join(osp.dirname(configPyFile), "../../data")

def data_path(fname):
    return osp.join(data_dir, fname)

gencodeBb = "gencodeV41.bb"

hg38gd = GenomeData("hg38",
                    data_path("hg38.2bit"),
                    assembly_report=data_path("GCF_000001405.40_GRCh38.p14_assembly_report.txt"))
hg38gd.add_track(osp.basename(gencodeBb),
                 data_path(gencodeBb),
                 f"https://hgdownload.soe.ucsc.edu/gbdb/hg38/gencode/${gencodeBb}")
hg38gd.add_track("WTC11_consolidated",
                 data_path("WTC11_consolidated.bigBed"),
                 "http://conesalab.org/LRGASP/LRGASP_hub/hg38/Human_samples/WTC11_consolidated.bigBed")

# no ispcr
hg38config = GenomeConfig(hg38gd)

primer3_config = Primer3Config()
primer3_config.PRIMER_NUM_RETURN = 12
primer3_config.thermo_args = Primer3ThermoArgs()
primer3_config.thermo_filters = DefaultThermoFitlers()

config = PrimersJuJuConfig(primer3_config=primer3_config)
config.add_genome(hg38config)
