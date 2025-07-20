from enum import Enum, auto
from pathlib import Path
from typing import Any, Callable

from pandas.core.frame import DataFrame
from xtraee.config import debug, get_logger
from xtraee.irrep import Irrep
from xtraee.lazypattern import LP
from xtraee.trblock import TransitionBlock

DatasetType = dict[str, dict[str, Any]]
PathType = str | Path

input_patterns = dict(
    ee_singlets=LP(r"^EE_SINGLETS\s*=?\s*\[?([\d,\s]+)\]?\s*"),
    ee_triplets=LP(r"^EE_TRIPLETS\s*=?\s*\[?([\d,\s]+)\]?\s*"),
)


class Block(Enum):
    null = auto()
    input = auto()
    trprops = auto()
    ee = auto()
    mo = auto()
    tr = auto()
    irrep = auto()


class BaseParser:
    name = "base"

    def __init__(
        self,
        input_file: PathType,
        threshold: float = 0.0,
    ):
        self.threshold = threshold
        self.input_file = input_file
        self.parser: dict[Block, Callable[[str], None]] = {
            Block.null: lambda line: None,
            Block.input: self.process_input,
        }
        self.log = get_logger(self.name, debug)
        self.reset()

    def reset(self) -> None:
        self.block = Block.null
        self.ee_singlets: list[int] = []
        self.N_singlets = 0
        self.ee_triplets: list[int] = []
        self.N_triplets = 0
        self.irreps: dict[str, Irrep] = {}
        self.data: DatasetType = {}
        self.homo: int = 0

    def detect_block(self, line: str) -> None:
        pass

    def process_input(self, line: str) -> None:
        # only valid for QChem inputs
        if self.N_singlets == 0:
            if (m := input_patterns["ee_singlets"].match(line)) is not None:
                self.ee_singlets = list(map(int, m.group(1).strip().split(",")))
                self.N_singlets = len(self.ee_singlets)
        if self.N_triplets == 0:
            if (m := input_patterns["ee_triplets"].match(line)) is not None:
                self.ee_triplets = list(map(int, m.group(1).strip().split(",")))
                self.N_triplets = len(self.ee_triplets)

    @staticmethod
    def select_lowest_excitations(states: list[Irrep]) -> list[TransitionBlock]:
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

    def write_full(self, fname: PathType) -> None:
        with open(fname, "w") as f:
            for irrep in self.irreps.values():
                f.write(f"{irrep.name}\n")
                for tr in irrep.trblocks:
                    f.write(f"{tr}\n")

    def write_dataset(self, fname: PathType) -> None:
        if len(self.data) == 0:
            self.data = data = self.create_dataset()
        with open(fname, "w") as f:
            f.write("---- Happy family ----\n")
            for s, d in data.items():
                f.write(f"------------{s}------------\n")
                f.write(f"{d['singlet']}\n")
                f.write("Matched triplets: " + "\n")  # noqa
                for tr, triplet in d["matched_triplets"].items():
                    f.write(f"{tr} <==> {',        '.join(triplet)}\n")
            f.write("--- End of the happy family :) ---")

    def create_dataset(self) -> DatasetType:
        for irr in self.irreps.values():
            irr.sort()
        singlets = []
        triplets = []
        for irr in self.irreps.values():
            if irr.multi == "singlet":
                singlets.append(irr)
            elif irr.multi == "triplet":
                triplets.append(irr)
            else:
                raise ValueError(f"Unknown excitation type {irr.multi}")
        lowest_singlets = BaseParser.select_lowest_excitations(singlets)
        data = {}
        for i, singlet in enumerate(lowest_singlets):
            matched_triplets = {}
            for idx, s_max in enumerate(singlet.transitions):
                if abs(s_max.amplitude) > self.threshold:
                    key = str(s_max)
                    matched_triplets[key] = []
                    for triplet in triplets:
                        for trblock in triplet.trblocks:
                            for tr in trblock.transitions:
                                if s_max.is_equal(tr):
                                    matched_triplets[key].append(
                                        f"{trblock.id_number} {trblock.irrep} {tr.amplitude:.4f} "
                                    )
                                    data[f"S{i + 1}"] = {
                                        "idx": idx,
                                        "singlet": singlet,
                                        "matched_triplets": matched_triplets,
                                    }
                                    break
        return data

    def compare_all(self, method: str) -> dict[str, DataFrame]:
        irreps_singlets = [
            irrep for irrep in self.irreps.values() if irrep.multi == "singlet"
        ]
        irreps_triplets = [
            irrep for irrep in self.irreps.values() if irrep.multi == "triplet"
        ]
        scores = {}
        for irrep_singlet in irreps_singlets:
            for irrep_triplet in irreps_triplets:
                scores[irrep_singlet.name + "_" + irrep_triplet.name] = (
                    irrep_singlet.compare(irrep_triplet, method)
                )
        return scores

    def match2std(self, my_irrep: Irrep, other_block: TransitionBlock):
        majors = []
        # match all of the states in the reference to current
        for i, trblock in enumerate(my_irrep.trblocks):
            if trblock.is_equal_std(other_block):
                majors.append(trblock)
        return majors

    def _find_equivalent_irrep(self, other_irrep: Irrep) -> Irrep:
        # first try respecting the name of the symmetry
        for irrep in self.irreps.values():
            if irrep.multi == other_irrep.multi and irrep.name == other_irrep.name:
                return irrep
        # if not found, try by multi
        for irrep in self.irreps.values():
            if irrep.multi == other_irrep.multi:
                return irrep
        # unreachable!
        raise ValueError(
            f"No equivalent irrep found for {other_irrep.multi} {other_irrep.name}"
        )

    def write_vs_std(self, path: str, o_irreps: DatasetType) -> None:
        data = {}
        for other_key, other_irrep in o_irreps.items():
            my_irrep = self._find_equivalent_irrep(other_irrep)
            data[other_key] = []
            for other_trblock in other_irrep.trblocks:
                result = {
                    "ref": other_trblock,
                    self.name: self.match2std(my_irrep, other_trblock),
                }
                data[other_key].append(result)
        with open(path, "w") as f:
            f.write("---- Sad family ----\n")
            for k, v in data.items():
                f.write(f"{k}\n")
                for item in v:
                    f.write(
                        f"ref: {item['ref'].multi}-{item['ref'].id_number}/{item['ref'].irrep} <=>"
                    )
                    line = ",".join(
                        map(
                            lambda block: f"{self.name}: {block.multi}-{block.id_number}/{block.irrep}",
                            item[self.name],
                        )
                    )
                    f.write(line + "\n")
            f.write("--- End of sad family :( ---")

    def compare_std(
        self, o_irreps: dict[str, Irrep], method: str
    ) -> dict[str, DataFrame]:
        irreps_singlets = [
            irrep for irrep in self.irreps.values() if irrep.multi == "singlet"
        ]
        irreps_triplets = [
            irrep for irrep in self.irreps.values() if irrep.multi == "triplet"
        ]
        oirr_singlets = [
            irrep for irrep in o_irreps.values() if irrep.multi == "singlet"
        ]
        oirr_triplets = [
            irrep for irrep in o_irreps.values() if irrep.multi == "triplet"
        ]
        scores = {}
        for irrep_singlet in irreps_singlets:
            for oirr_singlet in oirr_singlets:
                scores[irrep_singlet.identifier + "_vs_" + oirr_singlet.identifier] = (
                    irrep_singlet.compare(oirr_singlet, method)
                )
        for irrep_triplet in irreps_triplets:
            for oirr_triplet in oirr_triplets:
                scores[irrep_triplet.identifier + "_vs_" + oirr_triplet.identifier] = (
                    irrep_triplet.compare(oirr_triplet, method)
                )
        return scores
