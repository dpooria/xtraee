from typing import Dict, List

import pandas as pd

from xtraee.trblock import TransitionBlock


class Irrep:
    def __init__(self, name: str, ee_type: str, n_states: int, parent: str = ""):
        self.name = name
        self.n_states = n_states
        self.trblocks: List[TransitionBlock] = []
        self.trblocks_dict: Dict[str, TransitionBlock] = {}
        self.ee_type = ee_type
        self.sorted = False
        self.updated = True
        self.parent = parent

    @property
    def identifier(self) -> str:
        if self.parent:
            return f"{self.parent}-{self.name}.{self.ee_type}"
        else:
            return f"{self.name}.{self.ee_type}"

    def sort(self) -> None:
        self.sorted = True
        self.trblocks.sort(key=lambda t: t.excitation_energy)
        self.update_transitions()

    def scatter_attr(self, attr_name, attr_value):
        for tr in self.trblocks:
            setattr(tr, attr_name, attr_value)

    def update_transitions(self) -> None:
        if self.updated:
            return
        self.updated = True
        for t in self.trblocks:
            self.trblocks_dict[t.identifier] = t

    def append(self, transition: TransitionBlock) -> None:
        self.updated = False
        self.trblocks.append(transition)

    def __len__(self) -> int:
        return len(self.trblocks)

    def __getitem__(self, key) -> TransitionBlock:
        return self.trblocks[key]

    def compare(self, other, method: str) -> pd.DataFrame:
        scores = []
        rows = []
        for tr in self.trblocks:
            rows.append(tr.identifier)
            score_tr = []
            columns = []
            for o_tr in other.trblocks:
                columns.append(o_tr.identifier)
                score, *_ = tr.compare(o_tr, method)
                score_tr.append(score)
            scores.append(score_tr)
        return pd.DataFrame(scores, columns=columns, index=rows)
