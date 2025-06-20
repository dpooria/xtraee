import re
from enum import Enum
from typing import Dict, Optional

from pandas.core.frame import DataFrame
from xtraee.irrep import Irrep
from xtraee.parser.base import BaseParser, DatasetType, PathType
from xtraee.trblock import CISTransitionBlock, EOMEETransitionBlock, TransitionBlock

start_indicators = {"cisee": "CIS Excitation Energies", "mo": "Orbital Energies (a.u.)"}

trblock_begin_pattern = re.compile(
    r"^Excited\s+state\s+(\d+)\s*:\s*excitation\s+energy\s+\(eV\)\s*=\s*([-+]?\d+.\d+)\s*$"
)


class Block(Enum):
    null = 0
    ee = 1
    mo = 2


class QCISParser(BaseParser):
    name = "CIS"

    def __init__(
        self,
        input_file: PathType,
        threshold: float = 0.0,
    ):
        super().__init__(input_file, threshold)
        self.block = Block.null
        self.parser = {
            Block.null: lambda line: None,
            Block.ee: self.process_trblock,
            Block.mo: lambda line: None,
        }
        self.vsccsd: DatasetType = {}
        self._inside_eomee = False
        self._current_transition: Optional[TransitionBlock] = None
        self.irrep_singlets = Irrep("singlet", "singlet", 0)
        self.irrep_triplets = Irrep("triplet", "triplet", 0)
        self.homo = 0

    def detect_block(self, line: str) -> None:
        if start_indicators["cisee"] in line:
            self.block = Block.ee
        elif start_indicators["mo"] in line:
            self.process_irreps()
            self.block = Block.mo

    def process_trblock(self, line: str) -> None:
        if (m := trblock_begin_pattern.match(line)) is not None:
            if self._current_transition is not None:
                self._current_transition.sort()
                if self._current_transition.excitation == "singlet":
                    self.irrep_singlets.append(self._current_transition)
                elif self._current_transition.excitation == "triplet":
                    self.irrep_triplets.append(self._current_transition)
                else:
                    raise ValueError(f"Unknown excitation {self._current_transition}")
            self._current_transition = CISTransitionBlock(int(m.group(1)))
            self._current_transition.excitation_energy = float(m.group(2))
        elif self._current_transition is not None:
            self._current_transition.add_data(line)

    def process_irreps(self) -> None:
        self.irrep_singlets.sort()
        self.irrep_triplets.sort()
        self.irrep_singlets.n_states = len(self.irrep_singlets)
        self.irrep_triplets.n_states = len(self.irrep_triplets)
        singlet_trblocks = self.irrep_singlets.trblocks
        triplet_trblocks = self.irrep_triplets.trblocks
        self.homo = 0
        for tr in singlet_trblocks:
            self.homo = max(
                *[int(id[1]) for tr_ in tr.transitions for id in tr_.id_i], self.homo
            )
        for tr in triplet_trblocks:
            self.homo = max(
                *[int(id[1]) for tr_ in tr.transitions for id in tr_.id_i], self.homo
            )
        for tr in singlet_trblocks:
            tr.homo = self.homo
        for tr in triplet_trblocks:
            tr.homo = self.homo
        self.irreps_dict = {
            "singlet": self.irrep_singlets,
            "triplet": self.irrep_triplets,
        }

    def match2ccsd(self, key: str, ccsd_block: EOMEETransitionBlock):
        irrep = self.irreps_dict[key]
        majors = []
        # match all of the states in CCSD to CIS
        for i, trblock in enumerate(irrep.trblocks):
            if trblock.is_equal_eomee(ccsd_block):
                majors.append(trblock)
        return majors

    def compare_with_ccsd(
        self,
        ccsd_irrep_dict: DatasetType,
    ) -> DatasetType:
        data = {}
        for key_irrep_ccsd, irrep_ccsd in ccsd_irrep_dict.items():
            data[key_irrep_ccsd] = []
            for trblock_ccsd in irrep_ccsd.trblocks:
                result = {
                    "CCSD": trblock_ccsd,
                    "CIS": self.match2ccsd(trblock_ccsd.excitation, trblock_ccsd),
                }
                data[key_irrep_ccsd].append(result)
        return data

    def write_vsccsd(self, path: str, ccsd_irrep_dict: DatasetType) -> None:
        self.vsccsd = data = self.compare_with_ccsd(ccsd_irrep_dict)
        with open(path, "w") as f:
            f.write("---- Sad family ----\n")
            for k, v in data.items():
                f.write(f"{k}\n")
                for item in v:
                    f.write(
                        f"CCSD: {item['CCSD'].excitation}-{item['CCSD'].id_number}/{item['CCSD'].irrep} <=>"
                    )
                    line = ",".join(
                        map(
                            lambda block: f"CIS: {block.excitation}-{block.id_number}/{block.irrep}",
                            item["CIS"],
                        )
                    )
                    f.write(line + "\n")
            f.write("--- End of sad family :( ---")

    def compare_eomee(
        self, irreps_dict: Dict[str, Irrep], method: str
    ) -> Dict[str, DataFrame]:
        irreps_singlets = [
            irrep for irrep in self.irreps_dict.values() if irrep.ee_type == "singlet"
        ]
        irreps_triplets = [
            irrep for irrep in self.irreps_dict.values() if irrep.ee_type == "triplet"
        ]
        oirr_singlets = [
            irrep for irrep in irreps_dict.values() if irrep.ee_type == "singlet"
        ]
        oirr_triplets = [
            irrep for irrep in irreps_dict.values() if irrep.ee_type == "triplet"
        ]
        scores = {}
        for irrep_singlet in irreps_singlets:
            for oirr_singlet in oirr_singlets:
                scores["CIS_" + irrep_singlet.name + "_CCSD_" + oirr_singlet.name] = (
                    irrep_singlet.compare(oirr_singlet, method)
                )
        for irrep_triplet in irreps_triplets:
            for oirr_triplet in oirr_triplets:
                scores["CIS_" + irrep_triplet.name + "_CCSD_" + oirr_triplet.name] = (
                    irrep_triplet.compare(oirr_triplet, method)
                )
        return scores
