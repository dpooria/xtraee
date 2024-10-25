from typing import Dict, List

import pandas as pd

from xtraee.trblock import TransitionBlock


class Irrep:
    def __init__(
        self,
        name: str,
        ee_type: str,
        n_states: int,
    ):
        self.name = name
        self.n_states = n_states
        self.transitions: List[TransitionBlock] = []
        self.transitions_dict: Dict[str, TransitionBlock] = {}
        self.ee_type = ee_type

    def sort(self) -> None:
        self.transitions.sort(key=lambda t: t.excitation_energy)
        self.update_transitions()

    def update_transitions(self) -> None:
        for t in self.transitions:
            self.transitions_dict[t.identifier] = t

    def append(self, transition: TransitionBlock) -> None:
        self.transitions.append(transition)

    def __len__(self) -> int:
        return len(self.transitions)

    def __getitem__(self, key) -> TransitionBlock:
        return self.transitions[key]

    def compare(self, other, method: str) -> pd.DataFrame:
        scores = []
        rows = []
        for tr in self.transitions:
            rows.append(tr.identifier)
            score_tr = []
            columns = []
            for o_tr in other.transitions:
                columns.append(o_tr.identifier)
                score_tr.append(tr.compare(o_tr, method))
            scores.append(score_tr)
        return pd.DataFrame(scores, columns=columns, index=rows)
