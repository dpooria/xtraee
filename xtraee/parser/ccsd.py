import re
from typing import List, Dict, Optional, Any


class Transition:
    def __init__(
        self,
        amplitude: float,
        id_i: int | str,
        irr_i: str,
        s_i: str,
        id_f: int | str,
        irr_f: str,
        s_f: str,
        is_double: bool = False,
    ):
        self.amplitude = amplitude
        self.id_i = id_i
        self.irr_i = irr_i
        self.s_i = s_i
        self.id_f = id_f
        self.irr_f = irr_f
        self.s_f = s_f
        self.is_double = is_double

    def __and__(self, other) -> bool:
        if self.is_double:
            return (self.irr_i == other.irr_i) and (self.irr_f == other.irr_f)
        else:
            return (
                (self.id_i == other.id_i)
                and (self.irr_i == other.irr_i)
                and (self.id_f == other.id_f)
                and (self.irr_f == other.irr_f)
            )

    @classmethod
    def from_match(cls, match: re.Match):
        return cls(
            float(match.group(1)),
            int(match.group(2)),
            match.group(3),
            match.group(4),
            int(match.group(5)),
            match.group(6),
            match.group(7),
        )

    @classmethod
    def from_str(cls, line: str):
        s = line.split("->")
        s[1] = s[1].strip()
        if len(s) != 2:
            raise ValueError(f"cannot match {line}")
        # match with the float
        if (m := re.match(r"^\s*([-+]?\d+\.\d+)\s*.*$", s[0])) is not None:
            amp_str = m.group(1)
            amplitude = float(amp_str)
            # remove the match part
            idx = s[0].index(amp_str) + len(amp_str)
            s[0] = s[0][idx:].strip()
        else:
            print(f"Cannot match amplitude {line}")
            amplitude = 0.0
        # match with the transition
        return cls(
            amplitude,
            irr_i=s[0],
            irr_f=s[1],
            id_i="",
            id_f="",
            s_i="",
            s_f="",
            is_double=True,
        )

    def __repr__(self) -> str:
        return (
            f"{self.amplitude:.4f}\t{self.id_i} ({self.irr_i}) {self.s_i} -> "
            f"{self.id_f} ({self.irr_f}) {self.s_f}\n"
        )


class EOMEETransitionBlock:
    EE_PATTERN = re.compile(r"^.*Excitation energy\s*=\s*([-+]?\d*\.?\d+)\s*eV\.\s*$")
    R_PATTERN = re.compile(
        r"^.*R0\^2\s*=\s*(\d*.\d+)\s*R1\^2\s*=\s*([-+]?\d*\.?\d+)\s*R2\^2\s*=\s*([-+]?\d*\.?\d+).*$"
    )
    TRANSITION_PATTERN = re.compile(
        r"^\s*([-+]?\d+\.\d+)\s+(\d+)\s+\((.+)\)\s+(\w+)\s+->\s+(\d+)\s+\((.+)\)\s+(\w+)\s*$"
    )
    END_TRANSITIONBLOCK = "Summary of significant orbitals:"

    def __init__(self, id_number: int, irrep: str, excitation: str):
        self.transitions: List[Transition] = []
        self.id_number = id_number
        self.irrep = irrep
        self.excitation = excitation
        self.completed = False
        self.excitation_energy = 0.0
        self.R0 = 0.0
        self.R1 = 0.0
        self.R2 = 0.0
        self.identifier = f"{self.excitation}-{self.id_number}/{self.irrep}"
        self.oscillator_strength = 0.0
        self.omega = 0.0

    def add_data(self, line: str) -> bool:
        if self.END_TRANSITIONBLOCK in line:
            self.completed = True
            return True
        # if (m := self.TRANSITION_PATTERN.match(line)) is not None:
        #     self.transitions.append(Transition.from_match(m))
        # TODO: all the excitations are now double! fix it
        if "->" in line:
            self.transitions.append(Transition.from_str(line))
        elif (m := self.EE_PATTERN.match(line)) is not None:
            self.excitation_energy = float(m.group(1))
        elif (m := self.R_PATTERN.match(line)) is not None:
            self.R0 = float(m.group(1))
            self.R1 = float(m.group(2))
            self.R2 = float(m.group(3))
        return False

    def sort(self) -> None:
        self.transitions.sort(key=lambda t: abs(t.amplitude), reverse=True)

    def __repr__(self) -> str:
        return (
            f"EOMEE transition {self.excitation} {self.id_number}/{self.irrep}\n"  # noqa
            f"EE: {self.excitation_energy:.4f} eV.\n"
            f"R0^2: {self.R0:.4f} R1^2: {self.R1:.4f} R2^2: {self.R2:.4f}\n"  # noqa
            "Amplitude Transitions between orbitals\n"
            f'{"".join(map(str, self.transitions))}\n'
            f"Oscillator strength (a.u.): {self.oscillator_strength:.6f},"
            f"omega (Mulliken): {self.omega:.4f}\n"
        )

    # def __str__(self) -> str:
    #     return (
    #         f'EOMEE transition {self.excitation} {
    #             self.id_number}/{self.irrep}\n'
    #         f'EE: {self.excitation_energy:.4f} eV.\n'
    #         f'R0^2: {self.R0:.4f} R1^2: {
    #                 self.R1:.4f} R2^2: {self.R2:.4f}\n'
    #         f'Oscillator strength (a.u.): {self.oscillator_strength:.6f},'
    #         f'omega (Mulliken): {self.omega:.4f}\n')

    def __getitem__(self, key) -> Transition:
        return self.transitions[key]


class Irrep:
    def __init__(
        self,
        name: str,
        ee_type: str,
        n_states: int,
    ):
        self.name = name
        self.n_states = n_states
        self.eomee_transitions: List[EOMEETransitionBlock] = []
        self.transitions: Dict[str, EOMEETransitionBlock] = {}
        self.ee_type = ee_type

    def sort(self) -> None:
        self.eomee_transitions.sort(key=lambda t: t.excitation_energy)
        self.update_transitions()

    def update_transitions(self) -> None:
        for t in self.eomee_transitions:
            self.transitions[t.identifier] = t

    def append(self, transition: EOMEETransitionBlock) -> None:
        self.eomee_transitions.append(transition)

    def __len__(self) -> int:
        return len(self.eomee_transitions)

    def __getitem__(self, key) -> EOMEETransitionBlock:
        return self.eomee_transitions[key]


class CCSDParser:
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
    OMEGA_PATTERN = re.compile(r"^\s*omega\s+=\s+([-+]?\d+\.\d+)\s*$")

    IRREPSOLV_BEGIN = "Solving for EOMEE-CCSD"
    INPUT_BLOCK_BEGIN = "$rem"
    INPUT_BLOCK_END = "$end"
    LAMBDA_BLOCK_BEGIN = "CCSD Lambda converged."
    LAMBDA_BLOCK_END = "Start computing the transition properties"
    TRPROP_BLOCK_END = "All requested transition properties have been computed."

    NULL_BLOCK = ""
    INPUT_BLOCK = "input"
    LAMBDA_BLOCK = "lambda"
    TRPROPS_BLOCK = "transition_prop"

    def __init__(self, input_file: str, first_kid: str, happy_family: str):
        self.input_file = input_file
        self.first_kid = first_kid
        self.happy_family = happy_family
        self.parser = {
            CCSDParser.NULL_BLOCK: lambda line: None,
            CCSDParser.INPUT_BLOCK: self.process_input_block,
            CCSDParser.LAMBDA_BLOCK: self.process_lambda_block,
            CCSDParser.TRPROPS_BLOCK: self.process_trprops_block,
        }
        self.block = CCSDParser.NULL_BLOCK
        self.ee_singlets: List[int] = []
        self.N_singlets = 0
        self.ee_triplets: List[int] = []
        self.N_triplets = 0
        self.inside_eomee = False
        self.irreps: List[Irrep] = []
        self.irreps_dict: Dict[str, Irrep] = {}
        self.current_transition: Optional[EOMEETransitionBlock] = None
        self.current_irrep: Optional[Irrep] = None
        self.singlet_irrep_counter = 0
        self.triplet_irrep_counter = 0
        self.current_excitation: str = ""
        self.current_trprop: str = ""

    def process_file(self) -> None:
        with open(self.input_file, "r") as f:
            for line in f:
                self.parse_line(line)

    def parse_line(self, line: str) -> None:
        self.detect_block(line)
        self.parser[self.block](line)

    @staticmethod
    def select_lowest_excitations(states: List[Irrep]) -> List[EOMEETransitionBlock]:
        n_states = [irr.n_states for irr in states]
        lowest_excitations = []
        for i in range(min(n_states)):
            e_block = min([irr[i] for irr in states], key=lambda t: t.excitation_energy)
            lowest_excitations.append(e_block)
        return lowest_excitations

    # Transition | EOMEETransitionBlock | List[EOMEETransitionBlock]
    def create_happy_family(self) -> Dict[str, Dict[str, Any]]:
        for irr in self.irreps:
            irr.sort()
        singlets = [irr for irr in self.irreps if irr.ee_type == "singlet"]
        triplets = [irr for irr in self.irreps if irr.ee_type == "triplet"]
        lowest_singlets = CCSDParser.select_lowest_excitations(singlets)
        data = {}
        for i, singlet in enumerate(lowest_singlets):
            matched_triplets = []
            for idx, s_max in enumerate(singlet.transitions):
                for triplet in triplets:
                    for eomee_tr in triplet.eomee_transitions:
                        for tr in eomee_tr.transitions:
                            if s_max & tr:
                                matched_triplets.append(eomee_tr)
                                break
                if len(matched_triplets) > 0:
                    break
            data[f"S{i+1}"] = {
                "idx": idx,
                "s_max": s_max,
                "singlet": singlet,
                "matched_triplets": matched_triplets,
            }
        return data

    def write_first_kid(self) -> None:
        with open(self.first_kid, "w") as f:
            for irrep in self.irreps:
                f.write(f"{irrep.name}\n")
                for tr in irrep.eomee_transitions:
                    f.write(f"{tr}\n")

    def write_happy_family(self) -> None:
        data = self.create_happy_family()
        with open(self.happy_family, "w") as f:
            f.write("Lowest excitation transitions\n")
            f.write("---- Happy family -----\n")
            for s, d in data.items():
                f.write(f"------------{s}------------\n")
                f.write(f'{d["singlet"]}\n')
                f.write(
                    f'Matched triplets (with {d["idx"]}: {d["s_max"]} ' + "\n"
                )  # noqa
                for triplet in d["matched_triplets"]:
                    f.write(f"{triplet}\n")
            f.write("--- End of happy family :) -----\n")

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
                self.irreps.append(self.current_irrep)
                self.irreps_dict[f"{ee_type}-0/{m.group(1)}"] = (
                    self.current_irrep
                )  # noqa
                self.current_excitation = ee_type
        elif self.inside_eomee and self.current_transition is not None:
            if self.current_transition.add_data(line):
                self.current_transition.sort()
                self.inside_eomee = False
                if self.current_irrep is not None:
                    self.current_irrep.append(self.current_transition)
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

    def process_trprops_block(self, line: str) -> None:
        if (m := self.TRPROP_PATTERN.match(line)) is not None:
            ee_type, id_number, irrep = m.group(1), m.group(2), m.group(3)
            self.current_trprop = f"{ee_type}-{id_number}/{irrep}"
            self.current_irrep = self.irreps_dict[f"{ee_type}-0/{irrep}"]
            self.current_irrep.update_transitions()
            self.current_transition = self.current_irrep.transitions[
                self.current_trprop
            ]
        elif self.current_trprop != "" and self.current_transition is not None:
            if (m := self.OSCILLATOR_PATTERN.match(line)) is not None:
                self.current_transition.oscillator_strength = float(m.group(1))
            elif (m := self.OMEGA_PATTERN.match(line)) is not None:
                self.current_transition.omega = float(m.group(1))
                self.current_trprop = ""
