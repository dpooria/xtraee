from xtraee.irrep import Irrep
from xtraee.lazypattern import LP
from xtraee.parser.base import BaseParser, Block, PathType
from xtraee.trblock import CCSDTransitionBlock

meta_patterns = dict(
    irrepsolv=LP(
        r"^\s*Solving\s+for\s+EOMEE-(CCSD|CC2)\s+(.+)\s+(singlet|triplet)\s+states\.\s*$"
    ),
    eomee=LP(r"^\s*EOMEE\s+transition\s+(\d+)/(.+)\s*$"),
    trprop=LP(
        r"^\s*State\s+B:\s+eomee_ccsd/rhfref/(singlet|triplet)s:\s+(\d+)/(.+)\s*$"
    ),
    eeprop=LP(
        r"^\s*Excited state properties for\s+EOMEE-(CCSD|CC2) transition\s+(\d+)/(.+)\s*$"
    ),
)
prop_patterns = dict(
    oscillator_strength=LP(
        r"^\s*Oscillator strength \(a\.u\.\):\s+([-+]?\d+\.\d+)\s*$"
    ),
    gamma=LP(r"^\s*\|\|gamma\^AB\|\|\*\|\|gamma\^BA\|\|:\s*([0-9]+\.[0-9]+)\s*$"),
    omega=LP(r"^\s*omega\s+=\s+([-+]?\d+\.\d+)\s*$"),
    alphabeta=LP(r"^\s*2\<alpha\|beta\>\s+=\s+([-+]?\d+\.\d+)\s*$"),
    loc=LP(r"^\s*LOC\s+=\s+([-+]?\d+\.\d+)\s*$"),
    phe=LP(r"^\s*\<Phe\>\s+=\s+([-+]?\d+\.\d+)\s*$"),
    rhre=LP(r"^\s*\|<r_e - r_h>\|\s*\[Ang\]:\s*(\d+\.\d+)$"),
    corr_coef=LP(r"^\s*Correlation coefficient:\s*([-+]?[0-9]+\.[0-9]+)\s*$"),
)

start_indicators = dict(
    irrepsolv="Solving for",
    input="$rem",
    eeprop="Excited state properties for",
)
end_indicators = dict(
    input="$end",
    irrepsolv="Start computing the transition properties",
    trprop="All requested transition properties have been computed.",
)


class QCCSDParser(BaseParser):
    name = "EOM-CCSD"

    def __init__(
        self,
        input_file: PathType,
        threshold: float = 0.0,
    ):
        super().__init__(input_file, threshold)
        self.parser.update(
            {Block.irrep: self.process_irrepsolv, Block.trprops: self.process_trprops}
        )

    def reset(self):
        super().reset()
        self._inside_eomee = False
        self._inside_eeprop = False
        self._current_trblock: None | CCSDTransitionBlock = None
        self._current_irrep: None | Irrep = None
        self._current_multi: str = ""
        self._current_trprop: str = ""
        self._singlet_irrep_counter = 0
        self._triplet_irrep_counter = 0

    def detect_block(self, line: str) -> None:
        block = self.block
        if block == Block.null:
            if start_indicators["input"] in line:
                block = Block.input
            elif start_indicators["irrepsolv"] in line:
                block = Block.irrep
        elif block == Block.input and end_indicators["input"] in line:
            block = Block.null
        elif block == Block.irrep and end_indicators["irrepsolv"] in line:
            block = Block.trprops
        elif block == Block.null and end_indicators["trprop"] in line:
            block = Block.null
        self.block = block

    def process_irrepsolv(self, line: str) -> None:
        if start_indicators["irrepsolv"] in line:
            if (m := meta_patterns["irrepsolv"].match(line)) is not None:
                multi = m.group(3)
                if multi == "singlet":
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
                    m.group(2), multi, n_states, parent=self.name
                )
                self.irreps[f"{multi}-{m.group(2)}"] = self._current_irrep
                self._current_multi = multi
        elif self._inside_eomee and self._current_trblock is not None:
            if self._current_trblock.add_data(line):
                self._current_trblock.sort()
                self._inside_eomee = False
                if self._current_irrep is not None:
                    self._current_irrep.append(self._current_trblock)
        elif self._inside_eeprop:
            if not self._current_trblock.completed_extras:
                self._current_trblock.add_data(line)
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
                                irrep, self._current_irrep.name
                            )
                        )
                # else:
                #     raise ValueError("No current transition block")
                self._current_trblock = CCSDTransitionBlock(
                    int(m.group(1)), irrep, self._current_multi, name=self.name
                )
            elif start_indicators["eeprop"] in line:
                if (m := meta_patterns["eeprop"].match(line)) is not None:
                    self._inside_eeprop = True
                    irrep = m.group(3)
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
                    multi, id_number, irrep = (
                        self._current_irrep.multi,
                        m.group(2),
                        m.group(3),
                    )
                    self._current_trprop = f"{multi}-{id_number}/{irrep}"
                    self._current_trblock = self._current_irrep.trblocks_dict[
                        self._current_trprop
                    ]

    def process_trprops(self, line: str) -> None:
        if (m := meta_patterns["trprop"].match(line)) is not None:
            multi, id_number, irrep = m.group(1), m.group(2), m.group(3)
            self._current_trprop = f"{multi}-{id_number}/{irrep}"
            self._current_irrep = self.irreps[f"{multi}-{irrep}"]
            self._current_irrep.update_transitions()
            self._current_trblock = self._current_irrep.trblocks_dict[
                self._current_trprop
            ]
        elif self._current_trprop != "" and self._current_trblock is not None:
            for name, pattern in prop_patterns.items():
                if (m := pattern.match(line)) is not None:
                    value = float(m.group(1))
                    if name == "oscillator_strength":
                        self._current_trblock.oscillator_strength = value
                        continue
                    self._current_trblock.meta_data[name] = value
                    if name == "corr_coef":
                        self._current_trprop = ""  # end the current trprop state
                    break

    def gather_descriptors(self, id_prefix: str = ""):
        data = []
        for irr in self.irreps.values():
            irr.sort()
            for trblock in irr.trblocks:
                meta = trblock.meta_data.copy()
                meta.pop("R0")
                meta.pop("R1")
                meta["identifier"] = id_prefix + trblock.identifier
                meta["|r_e-r_h|"] = meta.pop("rhre")
                frontier_no = meta.pop("frontier_no")
                meta["frontier_no_1"] = frontier_no[0]
                meta["frontier_no_2"] = frontier_no[1]
                data.append(meta)

        return data


class QCC2Parser(QCCSDParser):
    name = "CC2"
