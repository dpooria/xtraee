import pandas as pd

from xtraee.trblock import TransitionBlock


class Irrep:
    def __init__(self, name: str, multi: str, n_states: int, parent: str = ""):
        self.name = name
        self.n_states = n_states
        self.trblocks: list[TransitionBlock] = []
        self.trblocks_dict: dict[str, TransitionBlock] = {}
        self.multi = multi
        self.sorted = False
        self.updated = True
        self.parent = parent

    @property
    def identifier(self) -> str:
        if self.parent:
            return f"{self.parent}-{self.name}.{self.multi}"
        else:
            return f"{self.name}.{self.multi}"

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

    def compare(self, other: "Irrep", method: str) -> pd.DataFrame:
        """Score every state of this irrep against every state of ``other``.

        Returns a state-by-state matrix.  When either side holds no parsed
        states the matrix is empty but still labelled, so ``compareall`` can
        write it instead of dying on an irrep that produced no states
        (previously ``UnboundLocalError: columns``).
        """
        rows = [tr.identifier for tr in self.trblocks]
        columns = [o_tr.identifier for o_tr in other.trblocks]
        if not rows or not columns:
            return pd.DataFrame(index=rows, columns=columns, dtype=float)
        scores = []
        for tr in self.trblocks:
            scores.append([tr.compare(o_tr, method)[0] for o_tr in other.trblocks])
        return pd.DataFrame(scores, columns=columns, index=rows)
