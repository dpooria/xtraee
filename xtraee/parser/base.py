from pathlib import Path
from abc import abstractmethod
from typing import Any, Callable, Dict, List
from enum import Enum

from xtraee.irrep import Irrep
from xtraee.trblock import TransitionBlock


class Block(Enum):
    null = 0


class Parser:
    def __init__(
        self,
        input_file: str | Path,
        threshold: float = 0.0,
        first_kid: str | Path = "first_kid.txt",
        happy_family: str | Path = "happy_family.txt",
    ):
        self.threshold = threshold
        self.input_file = input_file
        self.first_kid = first_kid
        self.happy_family_path = happy_family
        self.parser: Dict[int, Callable[[str], None]] = {Block.null: lambda line: None}
        self.block = Block.null  # remember to change this in the inherited classes
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
