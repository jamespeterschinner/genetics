from genetics.parser import parse_text
from helpers import horizontal_index
from helpers import reverse
from helpers import vertical_index

__all__ = ['Pedigree']


class Pedigree(object):


    @property
    def rows(self):
        return len(self._rows)

    @property
    def columns(self):
        return len(self._columns)

    @property
    def shape(self):
        return self._shape

    @property
    def affected(self):
        row_map = reverse(self._rows)
        column_map = reverse(self._columns)
        return tuple((row_map[row], column_map[column]) for row, column in self._affected)

    @property
    def roots(self):
        return tuple(observation for observation in self._indexed_observations.values()
                     if observation.root)

    def __init__(self, text):
        parsed_input = parse_text(text)
        self._array = parsed_input.array
        self._shape = parsed_input.shape
        self._indexed_observations = parsed_input.indexed_observations
        self._affected = parsed_input.affected_observations
        self._rows = parsed_input.generation_mapping.row
        self._columns = parsed_input.generation_mapping.column

    def __getitem__(self, index):
        row, column = index
        return self._array[self._rows[row]][self._columns[column]]

    def __setitem__(self, key, value):
        row, column = key
        self._array[self._rows[row]][self._columns[column]] = value

    def __repr__(self):
        """Create a nice console display
        """

        # Top index
        column_index = horizontal_index(self.shape[1], self._columns,
                                        len(str(self.shape[0])) + 1)

        # LHS Index
        row_index = vertical_index(self.shape[0], self._rows)
        # Combine LHS index with array
        result = ''.join(row + ' '.join(str(i) for i in self._array[idx])
                         for idx, row in enumerate(row_index))

        return f'{column_index}' \
               f'{result}'

    @classmethod
    def from_file(cls, file):
        with open(file, 'r') as f:
            return cls(f.read())
