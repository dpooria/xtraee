import re
from enum import Enum
from typing import Optional

from xtraee.irrep import Irrep
from xtraee.trblock import CCSDTransitionBlock, TransitionBlock

from xtraee.parser.base import BaseParser, PathType

meta_patterns = dict(
    ee_singlets=re.compile(r"^EE_SINGLETS\s+\[(.*)\]\s*"),
    ee_triplets=re.compile(r"^EE_TRIPLETS\s+\[(.*)\]\s*"),
    irrepsolv=re.compile(
        r"^\s*Solving\s+for\s+EOMEE-CCSD\s+(.+)\s+(singlet|triplet)\s+states\.\s*$"
    ),
    eomee=re.compile(r"^\s*EOMEE\s+transition\s+(\d+)/(.+)\s*$"),
    trprop=re.compile(
        r"^\s*State\s+B:\s+eomee_ccsd/rhfref/(singlet|triplet)s:\s+(\d+)/(.+)\s*$"
    ),
    eeprop=re.compile(
        r"^\s*Excited state properties for\s+EOMEE-CCSD transition\s+(\d+)/(.+)\s*$"
    ),
)
prop_patterns = dict(
    oscillator_strength=re.compile(
        r"^\s*Oscillator strength \(a\.u\.\):\s+([-+]?\d+\.\d+)\s*$"
    ),
    gamma=re.compile(
        r"^\s*\|\|gamma\^AB\|\|\*\|\|gamma\^BA\|\|:\s*([0-9]+\.[0-9]+)\s*$"
    ),
    omega=re.compile(r"^\s*omega\s+=\s+([-+]?\d+\.\d+)\s*$"),
    alphabeta=re.compile(r"^\s*2\<alpha\|beta\>\s+=\s+([-+]?\d+\.\d+)\s*$"),
    loc=re.compile(r"^\s*LOC\s+=\s+([-+]?\d+\.\d+)\s*$"),
    phe=re.compile(r"^\s*\<Phe\>\s+=\s+([-+]?\d+\.\d+)\s*$"),
    rhre=re.compile(r"^\s*\|<r_e - r_h>\|\s*\[Ang\]:\s*(\d+\.\d+)$"),
    corr_coef=re.compile(r"^\s*Correlation coefficient:\s*([-+]?[0-9]+\.[0-9]+)\s*$"),
)

start_indicators = dict(
    irrepsolv="Solving for EOMEE-CCSD",
    input="$rem",
    lambdab="CCSD Lambda converged.",
    eeprop="Excited state properties for  EOMEE-CCSD transition",
)
end_indicators = dict(
    input="$end",
    lambdab="Start computing the transition properties",
    trprop="All requested transition properties have been computed.",
)


class Block(Enum):
    null = 0
    input = 1
    lambdab = 2
    trprops = 3


class QCCSDParser(BaseParser):
    name = "EOM-CCSD"

    def __init__(
        self,
        input_file: PathType,
        threshold: float = 0.0,
    ):
        super().__init__(input_file, threshold)
        self.block = Block.null
        self.parser = {
            Block.null: lambda line: None,
            Block.input: self.process_input_block,
            Block.lambdab: self.process_lambda_block,
            Block.trprops: self.process_trprops_block,
        }
        self._inside_eomee = False
        self._inside_eeprop = False
        self._current_transition: Optional[TransitionBlock] = None
        self._current_irrep: Optional[Irrep] = None
        self._current_excitation: str = ""
        self._current_trprop: str = ""
        self._singlet_irrep_counter = 0
        self._triplet_irrep_counter = 0

    def detect_block(self, line: str) -> None:
        block = self.block
        if block == Block.null:
            if start_indicators["input"] in line:
                block = Block.input
            elif start_indicators["lambdab"] in line:
                block = Block.lambdab
        elif block == Block.input and end_indicators["input"] in line:
            block = Block.null
        elif block == Block.lambdab and end_indicators["lambdab"] in line:
            block = Block.trprops
        elif block == Block.null and end_indicators["trprop"] in line:
            block = Block.null
        self.block = block

    def process_input_block(self, line: str) -> None:
        if self.N_singlets == 0:
            if (m := meta_patterns["ee_singlets"].match(line)) is not None:
                self.ee_singlets = list(map(int, m.group(1).strip().split(",")))
                self.N_singlets = len(self.ee_singlets)
        if self.N_triplets == 0:
            if (m := meta_patterns["ee_triplets"].match(line)) is not None:
                self.ee_triplets = list(map(int, m.group(1).strip().split(",")))
                self.N_triplets = len(self.ee_triplets)

    def process_lambda_block(self, line: str) -> None:
        if start_indicators["irrepsolv"] in line:
            if (m := meta_patterns["irrepsolv"].match(line)) is not None:
                ee_type = m.group(2)
                if ee_type == "singlet":
                    n_states = self.ee_singlets[self._singlet_irrep_counter]
                    self._singlet_irrep_counter += 1
                else:
                    n_states = self.ee_triplets[self._triplet_irrep_counter]
                    self._triplet_irrep_counter += 1
                if self._current_irrep is not None:
                    cc = self._current_irrep
                    assert len(cc) == cc.n_states, "Irrep states mismatch"
                    cc.sort()
                self._current_irrep = Irrep(
                    m.group(1),
                    ee_type,
                    n_states,
                    parent="CCSD"
                )
                self.irreps_dict[f"{ee_type}-0/{m.group(1)}"] = self._current_irrep
                self._current_excitation = ee_type
        elif self._inside_eomee and self._current_transition is not None:
            if self._current_transition.add_data(line):
                self._current_transition.sort()
                self._inside_eomee = False
                if self._current_irrep is not None:
                    self._current_irrep.append(self._current_transition)
        elif self._inside_eeprop:
            if not self._current_transition.completed_extras:
                self._current_transition.add_data(line)
            else:
                self._inside_eeprop = False
        else:
            if (m := meta_patterns["eomee"].match(line)) is not None:
                self._inside_eomee = True
                irrep = m.group(2)
                if self._current_irrep is not None:
                    if self._current_irrep.name != irrep:
                        raise ValueError(
                            "Transition block irrep mismatch {} != {}".format(
                                irrep, self._current_irrep
                            )
                        )
                else:
                    raise ValueError("No current transition block")
                self._current_transition = CCSDTransitionBlock(
                    int(m.group(1)), irrep, self._current_excitation
                )
            elif start_indicators["eeprop"] in line:
                if (m := meta_patterns["eeprop"].match(line)) is not None:
                    self._inside_eeprop = True
                    irrep = m.group(2)
                    if self._current_irrep is not None:
                        if self._current_irrep.name != irrep:
                            raise ValueError(
                                "Transition block irrep mismatch {} != {}".format(
                                    irrep, self._current_irrep
                                )
                            )
                    else:
                        raise ValueError("No current irreducible representation")
                    self._current_irrep.update_transitions()
                    ee_type, id_number, irrep = (
                        self._current_irrep.ee_type,
                        m.group(1),
                        m.group(2),
                    )
                    self._current_trprop = f"{ee_type}-{id_number}/{irrep}"
                    self._current_transition = self._current_irrep.trblocks_dict[
                        self._current_trprop
                    ]

    def process_trprops_block(self, line: str) -> None:
        if (m := meta_patterns["trprop"].match(line)) is not None:
            ee_type, id_number, irrep = m.group(1), m.group(2), m.group(3)
            self._current_trprop = f"{ee_type}-{id_number}/{irrep}"
            self._current_irrep = self.irreps_dict[f"{ee_type}-0/{irrep}"]
            self._current_irrep.update_transitions()
            self._current_transition = self._current_irrep.trblocks_dict[
                self._current_trprop
            ]
        elif self._current_trprop != "" and self._current_transition is not None:
            for name, pattern in prop_patterns.items():
                if (m := pattern.match(line)) is not None:
                    value = float(m.group(1))
                    setattr(self._current_transition, name, value)
                    if name == "corr_coef":
                        self._current_trprop = ""  # end the current trprop state
                    break

    def gather_descriptors(self, id_prefix: str = ""):
        data = []
        for irr in self.irreps_dict.values():
            irr.sort()
            for trblock in irr.trblocks:
                data.append(
                    {
                        "id": id_prefix + trblock.identifier,
                        "R2": trblock.R2,
                        "gamma": trblock.gamma,
                        "omega": trblock.omega,
                        "loc": trblock.loc,
                        "phe": trblock.phe,
                        "alphabeta": trblock.alphabeta,
                        "|r_e-r_h|": trblock.rhre,
                        "corr_coef": trblock.corr_coef,
                        "froniter_no_1": trblock.froniter_no[0],
                        "froniter_no_2": trblock.froniter_no[1],
                        "nu": trblock.nu,
                        "nl": trblock.nl,
                        "prno": trblock.prno,
                    }
                )
        return data
