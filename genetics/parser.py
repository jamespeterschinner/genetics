from collections import namedtuple

from core import VALID_OBSERVATIONS, VALID_INSTRUCTIONS
from helpers import reverse
from observation import Observation

__all__ = ['parse_text']

ParsedText = namedtuple('ParsedText', ['array', 'shape'])

Index = namedtuple('Index', ['generation', 'child'])

ReifiedResults = namedtuple('ReifiedResults', ['indexed_observations', 'affected_observations', 'commands'])

GenerationMapping = namedtuple('GenerationMapping', ['row', 'column'])

ParsedInput = namedtuple('ParsedInput',
                         ['array', 'indexed_observations', 'shape', 'affected_observations', 'generation_mapping'])


def _create_text_array(text):
    rows = text.split('\n')
    # Pad out each row to the max len row
    max_len = max(map(len, rows))

    array = [list(str(row).ljust(max_len, ' ')) for row in rows]

    shape = (len(array), max_len)

    return ParsedText(array, shape)


def _reify_observations(array, typ, valid_observations, valid_instructions):
    """Reify each element into an Observation

    Args:
        array: numpy array
        typ: the type to reify each valid observation into
        valid_observations: A set of valid observations. Eg 'f' or 'M'
        valid_instructions: A set of valid instructions. Eg, '-' or '|'

    Returns:"""
    indexed_observations = {}
    affected_observations = []
    commands = []

    for row_idx, row in enumerate(array):
        for column_idx, character in enumerate(row):

            if character in valid_instructions:
                commands.append((row_idx, column_idx, character))

            elif character in valid_observations:
                observation = typ(character)
                indexed_observations.update({(row_idx, column_idx): observation})
                array[row_idx][column_idx] = observation

                if observation.affected:
                    affected_observations.append((row_idx, column_idx))

    return ReifiedResults(indexed_observations, affected_observations, commands)


def _link_partners(array, skip, commands, link_symbol):
    """Link partners in the pedigree

    Args:
        array: numpy array containing the pedigree information
        skip: a set of array indexes to skip linking
        commands: A list of commands itemizing what to link
        link_symbol: The symbol to indicate a link

    Returns:
        skip
    """

    for row, column, instruction in filter(lambda x: x[2] == link_symbol, reversed(commands)):
        partner_column = column + 1
        array[row][column - 1].partner = array[row][partner_column]
        skip.add((row, partner_column))  # Stop partners from being added to
        # the children of the same mother/father
    return skip


def _link_children(array, shape, skip, commands, link_symbol, typ):
    """Link children to parents

    Args:
        array: numpy array
        skip: a set of array indexes to skip
        commands: A list of commands itemizing what to link
        link_symbol: the symbol to indicate what to link
        typ: the type of a valid observation

    Returns:
        skip
    """

    # Link children
    for row, column, instruction in filter(lambda x: x[2] == link_symbol, reversed(commands)):
        parent = array[row - 1][column]
        mother, father = parent.mother_father
        child_row = row + 1
        for child_idx in range(column, shape[1]):
            if (child_row, child_idx) not in skip or child_idx == column:
                skip.add((child_row, child_idx))
                observation = array[child_row][child_idx]
                if isinstance(observation, typ):
                    observation.mother = mother
                    observation.father = father
    return skip


def _create_generation_mapping(indexed_observations):
    """Create individual mappings for generational indexes
     that map row and column values to the array index

     Args:
         indexed_observations: A mapping of the array index to to the observations
    """
    generation_rows = dict(enumerate(sorted(set(map(lambda x: x[0], indexed_observations.keys())))))
    generation_columns = dict(enumerate(sorted(set(map(lambda x: x[1], indexed_observations.keys())))))
    return GenerationMapping(generation_rows, generation_columns)


def _add_generation_index(indexed_observations, generation_rows, generation_columns):
    """Add the generation index to each observation

    Args:
        indexed_observations: A mapping of the array index to to the observations
        generation_rows: A mapping from array index row to generation index row
        generation_columns: A mapping from array index column to generation row
    """
    # generation_rows = dict(reverse(generation_rows))
    # generation_columns = dict(reverse(generation_columns))
    generation_rows, generation_columns = reverse(generation_rows), reverse(generation_columns)
    for (row, column), observation in indexed_observations.items():
        observation.generation = (generation_rows[row], generation_columns[column])


def parse_text(text):
    array, shape = _create_text_array(text)

    # Used to indicate that subsequent operations on each observation
    # should be skipped
    skip = set()

    # Instantiate Observation instances form the array
    reified_result = _reify_observations(
        array=array,
        typ=Observation,
        valid_observations=VALID_OBSERVATIONS,
        valid_instructions=VALID_INSTRUCTIONS
    )

    # Create the partner relationships between observations
    skip = _link_partners(
        array=array,
        skip=skip,
        commands=reified_result.commands,
        link_symbol='-'
    )

    # Create the children relationships between observations
    _link_children(
        array=array,
        shape=shape,
        skip=skip,
        commands=reified_result.commands,
        link_symbol='|',
        typ=Observation,
    )

    # A generation index is designed to indexes that correspond to 'Third generation fourth child'
    generation_mapping = _create_generation_mapping(
        indexed_observations=reified_result.indexed_observations
    )

    # Add each observations generation to the observation so that the observation
    # knows where it sits in the pedigree
    _add_generation_index(
        indexed_observations=reified_result.indexed_observations,
        generation_rows=generation_mapping.row,
        generation_columns=generation_mapping.column
    )

    return ParsedInput(array=array,
                       shape=shape,
                       indexed_observations=reified_result.indexed_observations,

                       affected_observations=reified_result.affected_observations,
                       generation_mapping=generation_mapping)
