from xtraee.trblock import TransitionBlock, EOMEETransitionBlock
from typing import List, Dict


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
