import re

VERBOSE = re.VERBOSE


class LP:
    __slots__ = ("_text", "_flags", "_pat")

    def __init__(self, text, flags=0):
        self._text = text
        self._flags = flags
        self._pat = None

    def _compile(self):
        if self._pat is None:
            self._pat = re.compile(self._text, self._flags)

    def search(self, s, *args, **kwargs):
        self._compile()
        return self._pat.search(s)

    def match(self, s, *args, **kwargs):
        self._compile()
        return self._pat.match(s, *args, **kwargs)

    def findall(self, s, *args, **kwargs):
        self._compile()
        return self._pat.findall(s, *args, **kwargs)

    def finditer(self, s, *args, **kwargs):
        self._compile()
        return self._pat.finditer(s, *args, **kwargs)

    def sub(self, sb, line, *args, **kwargs):
        self._compile()
        return self._pat.sub(sb, line, *args, **kwargs)


flp = LP(r"[-+]?\d+\.\d+")
slp = LP(r"\s+")
