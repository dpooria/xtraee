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

    def search(self, s):
        self._compile()
        return self._pat.search(s)

    def match(self, s):
        self._compile()
        return self._pat.match(s)

    def findall(self, s):
        self._compile()
        return self._pat.findall(s)

    def finditer(self, s):
        self._compile()
        return self._pat.finditer(s)
