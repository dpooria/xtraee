import pathlib
import re
from abc import abstractmethod
from typing import Any, Callable, Dict, List, Optional

from xtraee.irrep import Irrep
from xtraee.trblock import CISTransitionBlock, EOMEETransitionBlock, TransitionBlock


class Parser:
    NULL_BLOCK = 0

    def __init__(
        self,
        input_file: str | pathlib.Path,
        threshold: float = 0.0,
        first_kid: str | pathlib.Path = "first_kid.txt",
        happy_family: str | pathlib.Path = "happy_family.txt",
    ):
        self.threshold = threshold
        self.input_file = input_file
        self.first_kid = first_kid
        self.happy_family_path = happy_family
        self.parser: Dict[int, Callable[[str], None]] = {
            self.NULL_BLOCK: lambda line: None
        }
        self.block = Parser.NULL_BLOCK
        self.ee_singlets: List[int] = []
        self.N_singlets = 0
        self.ee_triplets: List[int] = []
        self.N_triplets = 0
        self.irreps_dict: Dict[str, Irrep] = {}

    @abstractmethod
    def detect_block(self, line: str) -> None:
        pass

    @staticmethod
    def select_lowest_excitations(states: List[Irrep]) -> List[TransitionBlock]:
        n_states = [irr.n_states for irr in states]
        lowest_excitations = []
        for i in range(min(n_states)):
            e_block = min([irr[i] for irr in states], key=lambda t: t.excitation_energy)
            lowest_excitations.append(e_block)
        return lowest_excitations

    def process_file(self) -> None:
        with open(self.input_file, "r") as f:
            for line in f:
                self.parse_line(line.strip())

    def parse_line(self, line: str) -> None:
        self.detect_block(line)
        self.parser[self.block](line)

    def write_first_kid(self) -> None:
        with open(self.first_kid, "w") as f:
            for irrep in self.irreps_dict.values():
                f.write(f"{irrep.name}\n")
                for tr in irrep.trblocks:
                    f.write(f"{tr}\n")

    def write_happy_family(self) -> None:
        self.happy_family = data = self.create_happy_family()
        with open(self.happy_family_path, "w") as f:
            f.write("---- Happy family ----\n")
            for s, d in data.items():
                f.write(f"------------{s}------------\n")
                f.write(f'{d["singlet"]}\n')
                f.write("Matched triplets: " + "\n")  # noqa
                for tr, triplet in d["matched_triplets"].items():
                    f.write(f'{tr} <==> {",        ".join(triplet)}\n')
            f.write("--- End of the happy family :) ---")

    def create_happy_family(self) -> Dict[str, Dict[str, Any]]:
        for irr in self.irreps_dict.values():
            irr.sort()
        singlets = []
        triplets = []
        for irr in self.irreps_dict.values():
            if irr.ee_type == "singlet":
                singlets.append(irr)
            elif irr.ee_type == "triplet":
                triplets.append(irr)
            else:
                raise ValueError(f"Unknown excitation type {irr.ee_type}")
        lowest_singlets = Parser.select_lowest_excitations(singlets)
        data = {}
        for i, singlet in enumerate(lowest_singlets):
            matched_triplets = {}
            for idx, s_max in enumerate(singlet.transitions):
                if abs(s_max.amplitude) > self.threshold:
                    ss = str(s_max).strip()
                    matched_triplets[ss] = []
                    for triplet in triplets:
                        for trblock in triplet.trblocks:
                            for tr in trblock.transitions:
                                if s_max.is_equal(tr):
                                    matched_triplets[ss].append(
                                        f"{trblock.id_number} {trblock.irrep} {tr.amplitude:.4f} "
                                    )
                                    break
            data[f"S{i + 1}"] = {
                "idx": idx,
                "singlet": singlet,
                "matched_triplets": matched_triplets,
            }
        return data

    def compare_all(self, method: str):
        irreps_singlets = [
            irrep for irrep in self.irreps_dict.values() if irrep.ee_type == "singlet"
        ]
        irreps_triplets = [
            irrep for irrep in self.irreps_dict.values() if irrep.ee_type == "triplet"
        ]
        scores = {}
        for irrep_singlet in irreps_singlets:
            for irrep_triplet in irreps_triplets:
                scores[irrep_singlet.name + "_" + irrep_triplet.name] = (
                    irrep_singlet.compare(irrep_triplet, method)
                )
        return scores


class QCCSDParser(Parser):
    EE_SINGLETS_PATTERN = re.compile(r"^EE_SINGLETS\s+\[(.*)\]\s*")
    EE_TRIPLETS_PATTERN = re.compile(r"^EE_TRIPLETS\s+\[(.*)\]\s*")
    IRREPSOLV_PATTERN = re.compile(
        r"^\s*Solving\s+for\s+EOMEE-CCSD\s+(.+)\s+(singlet|triplet)\s+states\.\s*$"
    )
    EOMEE_PATTERN = re.compile(r"^\s*EOMEE\s+transition\s+(\d+)/(.+)\s*$")
    TRPROP_PATTERN = re.compile(
        r"^\s*State\s+B:\s+eomee_ccsd/rhfref/(singlet|triplet)s:\s+(\d+)/(.+)\s*$"
    )
    OSCILLATOR_PATTERN = re.compile(
        r"^\s*Oscillator strength \(a\.u\.\):\s+([-+]?\d+\.\d+)\s*$"
    )
    GAMMA_PATTERN = re.compile(
        r"^\s*\|\|gamma\^AB\|\|\*\|\|gamma\^BA\|\|:\s*([0-9]+\.[0-9]+)\s*$"
    )
    OMEGA_PATTERN = re.compile(r"^\s*omega\s+=\s+([-+]?\d+\.\d+)\s*$")
    ALPHA_BETA_PATTERN = re.compile(r"^\s*2\<alpha\|beta\>\s+=\s+([-+]?\d+\.\d+)\s*$")
    LOC_PATTERN = re.compile(r"^\s*LOC\s+=\s+([-+]?\d+\.\d+)\s*$")
    CORRC_PATTERN = re.compile(
        r"^\s*Correlation coefficient:\s*([-+]?[0-9]+\.[0-9]+)\s*$"
    )
    EEPROP_PATTERN = re.compile(
        r"^\s*Excited state properties for\s+EOMEE-CCSD transition\s+(\d+)/(.+)\s*$"
    )

    IRREPSOLV_BEGIN = "Solving for EOMEE-CCSD"
    INPUT_BLOCK_BEGIN = "$rem"
    INPUT_BLOCK_END = "$end"
    LAMBDA_BLOCK_BEGIN = "CCSD Lambda converged."
    LAMBDA_BLOCK_END = "Start computing the transition properties"
    TRPROP_BLOCK_END = "All requested transition properties have been computed."
    EEPROP_BEGIN = "Excited state properties for  EOMEE-CCSD transition"

    INPUT_BLOCK = 1
    LAMBDA_BLOCK = 2
    TRPROPS_BLOCK = 3

    def __init__(
        self,
        input_file: str | pathlib.Path,
        threshold: float = 0.0,
        first_kid: str | pathlib.Path = "first_kid.txt",
        happy_family: str | pathlib.Path = "happy_family.txt",
    ):
        super().__init__(input_file, threshold, first_kid, happy_family)
        self.parser = {
            self.NULL_BLOCK: lambda line: None,
            self.INPUT_BLOCK: self.process_input_block,
            self.LAMBDA_BLOCK: self.process_lambda_block,
            self.TRPROPS_BLOCK: self.process_trprops_block,
        }
        self.inside_eomee = False
        self.inside_eeprop = False
        self.current_transition: Optional[TransitionBlock] = None
        self.current_irrep: Optional[Irrep] = None
        self.singlet_irrep_counter = 0
        self.triplet_irrep_counter = 0
        self.current_excitation: str = ""
        self.current_trprop: str = ""

    def detect_block(self, line: str) -> None:
        if self.block == self.NULL_BLOCK:
            if self.INPUT_BLOCK_BEGIN in line:
                self.block = self.INPUT_BLOCK
            elif self.LAMBDA_BLOCK_BEGIN in line:
                self.block = self.LAMBDA_BLOCK
        elif self.block == self.INPUT_BLOCK and self.INPUT_BLOCK_END in line:
            self.block = self.NULL_BLOCK
        elif self.block == self.LAMBDA_BLOCK and self.LAMBDA_BLOCK_END in line:
            self.block = self.TRPROPS_BLOCK
        elif self.block == self.TRPROPS_BLOCK and self.TRPROP_BLOCK_END in line:
            self.block = self.NULL_BLOCK

    def process_input_block(self, line: str) -> None:
        if self.N_singlets == 0:
            if (m := self.EE_SINGLETS_PATTERN.match(line)) is not None:
                self.ee_singlets = list(map(int, m.group(1).strip().split(",")))
                self.N_singlets = len(self.ee_singlets)
        if self.N_triplets == 0:
            if (m := self.EE_TRIPLETS_PATTERN.match(line)) is not None:
                self.ee_triplets = list(map(int, m.group(1).strip().split(",")))
                self.N_triplets = len(self.ee_triplets)

    def process_lambda_block(self, line: str) -> None:
        if self.IRREPSOLV_BEGIN in line:
            if (m := self.IRREPSOLV_PATTERN.match(line)) is not None:
                ee_type = m.group(2)
                if ee_type == "singlet":
                    n_states = self.ee_singlets[self.singlet_irrep_counter]
                    self.singlet_irrep_counter += 1
                else:
                    n_states = self.ee_triplets[self.triplet_irrep_counter]
                    self.triplet_irrep_counter += 1
                if self.current_irrep is not None:
                    cc = self.current_irrep
                    assert len(cc) == cc.n_states, "Irrep states mismatch"
                    cc.sort()
                self.current_irrep = Irrep(
                    m.group(1),
                    ee_type,
                    n_states,
                )
                self.irreps_dict[f"{ee_type}-0/{m.group(1)}"] = self.current_irrep
                self.current_excitation = ee_type
        elif self.inside_eomee and self.current_transition is not None:
            if self.current_transition.add_data(line):
                self.current_transition.sort()
                self.inside_eomee = False
                if self.current_irrep is not None:
                    self.current_irrep.append(self.current_transition)
        elif self.inside_eeprop:
            if not self.current_transition.completed_extras:
                self.current_transition.add_data(line)
            else:
                self.inside_eeprop = False
        else:
            if (m := self.EOMEE_PATTERN.match(line)) is not None:
                self.inside_eomee = True
                irrep = m.group(2)
                if self.current_irrep is not None:
                    if self.current_irrep.name != irrep:
                        raise ValueError(
                            "Transition block irrep mismatch {} != {}".format(
                                irrep, self.current_irrep
                            )
                        )
                else:
                    raise ValueError("No current transition block")
                self.current_transition = EOMEETransitionBlock(
                    int(m.group(1)), irrep, self.current_excitation
                )
            elif self.EEPROP_BEGIN in line:
                if (m := self.EEPROP_PATTERN.match(line)) is not None:
                    self.inside_eeprop = True
                    irrep = m.group(2)
                    if self.current_irrep is not None:
                        if self.current_irrep.name != irrep:
                            raise ValueError(
                                "Transition block irrep mismatch {} != {}".format(
                                    irrep, self.current_irrep
                                )
                            )
                    else:
                        raise ValueError("No current irreducible representation")
                    self.current_irrep.update_transitions()
                    ee_type, id_number, irrep = (
                        self.current_irrep.ee_type,
                        m.group(1),
                        m.group(2),
                    )
                    self.current_trprop = f"{ee_type}-{id_number}/{irrep}"
                    self.current_transition = self.current_irrep.trblocks_dict[
                        self.current_trprop
                    ]

    def process_trprops_block(self, line: str) -> None:
        if (m := self.TRPROP_PATTERN.match(line)) is not None:
            ee_type, id_number, irrep = m.group(1), m.group(2), m.group(3)
            self.current_trprop = f"{ee_type}-{id_number}/{irrep}"
            self.current_irrep = self.irreps_dict[f"{ee_type}-0/{irrep}"]
            self.current_irrep.update_transitions()
            self.current_transition = self.current_irrep.trblocks_dict[
                self.current_trprop
            ]
        elif self.current_trprop != "" and self.current_transition is not None:
            if (m := self.OSCILLATOR_PATTERN.match(line)) is not None:
                self.current_transition.oscillator_strength = float(m.group(1))
            elif (m := self.GAMMA_PATTERN.match(line)) is not None:
                assert isinstance(self.current_transition, EOMEETransitionBlock)
                self.current_transition.gamma = float(m.group(1))
            elif (m := self.OMEGA_PATTERN.match(line)) is not None:
                assert isinstance(self.current_transition, EOMEETransitionBlock)
                self.current_transition.omega = float(m.group(1))
            elif (m := self.ALPHA_BETA_PATTERN.match(line)) is not None:
                assert isinstance(self.current_transition, EOMEETransitionBlock)
                self.current_transition.alphabeta = float(m.group(1))
            elif (m := self.LOC_PATTERN.match(line)) is not None:
                assert isinstance(self.current_transition, EOMEETransitionBlock)
                self.current_transition.loc = float(m.group(1))
            elif (m := self.CORRC_PATTERN.match(line)) is not None:
                self.current_transition.corr_coef = float(m.group(1))
                # change this if you are going to extract more data from trprop
                # end the current transition this is the last value we extract from trprop
                self.current_trprop = ""

    def gather_descriptors(self, extra_id=""):
        data = []
        for irr in self.irreps_dict.values():
            irr.sort()
            for trblock in irr.trblocks:
                data.append(
                    {
                        "id": extra_id + trblock.identifier,
                        "R2": trblock.R2,
                        "gamma": trblock.gamma,
                        "omega": trblock.omega,
                        "loc": trblock.loc,
                        "alphabeta": trblock.alphabeta,
                        "corr_coef": trblock.corr_coef,
                        "froniter_no_1": trblock.froniter_no[0],
                        "froniter_no_2": trblock.froniter_no[1],
                        "nu": trblock.nu,
                        "nl": trblock.nl,
                        "prno": trblock.prno,
                    }
                )
        return data


class QCISParser(Parser):
    CISEE_BEGIN = "CIS Excitation Energies"
    MO_BEGIN = "Orbital Energies (a.u.)"
    TRBLOCK_BEGIN_PATTERN = re.compile(
        r"^Excited\s+state\s+(\d+)\s*:\s*excitation\s+energy\s+\(eV\)\s*=\s*([-+]?\d+.\d+)\s*$"
    )
    EE_BLOCK = 1
    MO_BLOCK = 2

    def __init__(
        self,
        input_file: str | pathlib.Path,
        threshold: float = 0.0,
        first_kid: str | pathlib.Path = "first_kid.txt",
        happy_family: str | pathlib.Path = "happy_family.txt",
    ):
        super().__init__(input_file, threshold, first_kid, happy_family)
        self.parser = {
            self.NULL_BLOCK: lambda line: None,
            self.EE_BLOCK: self.process_trblock,
            self.MO_BLOCK: lambda line: None,
        }
        self.sad_family: Dict[str, Dict[str, Any]] = {}
        self.inside_eomee = False
        self.current_transition: Optional[TransitionBlock] = None
        self.irrep_singlets = Irrep("singlet", "singlet", 0)
        self.irrep_triplets = Irrep("triplet", "triplet", 0)
        self.homo = 0

    def detect_block(self, line: str) -> None:
        if self.CISEE_BEGIN in line:
            self.block = self.EE_BLOCK
        elif self.MO_BEGIN in line:
            self.process_irreps()
            self.block = self.MO_BLOCK

    def process_trblock(self, line: str) -> None:
        if (m := self.TRBLOCK_BEGIN_PATTERN.match(line)) is not None:
            if self.current_transition is not None:
                self.current_transition.sort()
                if self.current_transition.excitation == "singlet":
                    self.irrep_singlets.append(self.current_transition)
                elif self.current_transition.excitation == "triplet":
                    self.irrep_triplets.append(self.current_transition)
                else:
                    raise ValueError(f"Unknown excitation {self.current_transition}")
            self.current_transition = CISTransitionBlock(int(m.group(1)))
            self.current_transition.excitation_energy = float(m.group(2))
        elif self.current_transition is not None:
            self.current_transition.add_data(line)

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

    def create_sad_family(
        self,
        ccsd_irrep_dict: Dict[str, Dict[str, Any]],
    ) -> Dict[str, Dict[str, Any]]:
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

    def write_sad_family(self, ccsd_irrep_dict: Dict[str, Dict[str, Any]], path: str):
        self.sad_family = data = self.create_sad_family(ccsd_irrep_dict)
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

    def compare_eomee(self, irreps_dict: Dict[str, Irrep], method: str):
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
