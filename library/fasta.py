from library.record import *
from typing import TextIO, Iterator, List, Tuple, Callable


def split_file(file: TextIO) -> Iterator[List[str]]:
    """
    Returns iterator that yield records as lists of lines
    """
    # find the beginning of the first record
    line = " "
    while line[0] != '>':
        line = file.readline()

    # chunk contains the already read lines of the current record
    chunk = []
    # put the first line of the first record into chunk
    chunk.append(line.rstrip())

    for line in file:
        # skip the blank lines
        if line == "" or line.isspace():
            continue

        # yield the chunk if the new record has begun
        if line[0] == '>':
            yield chunk
            chunk = []

        # put the first line of the new record into chunk
        chunk.append(line.rstrip())

    # yield the last record
    yield chunk


class Fastafile:
    """ Class for standard FASTA files"""

    @staticmethod
    def read(file: TextIO) -> Tuple[List[str], Callable[[], Iterator[Record]]]:
        """FASTA reader method"""

        # FASTA always have the same fields
        fields = ['seqid', 'sequence']

        def record_generator() -> Iterator[Record]:
            for chunk in split_file(file):
                # 'seqid' is the first line without the initial character
                # 'sequence' is the concatenation of all the other lines
                yield Record(seqid=chunk[0][1:], sequence="".join(chunk[1:]))
        return fields, record_generator
