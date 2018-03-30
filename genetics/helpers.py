
def reverse(obj):
    try:
        return reversed(obj)
    except TypeError:
        if isinstance(obj, dict):
            return {v: k for k, v in obj.items()}


def horizontal_index(n, mapping, offset=0):
    """Create horizontal index with vertically aligned numbers.
    Aligns numbers with mapping positions.

    Args:
        n: number of indexes to create
        mapping: A mapping of left->right index to equivalent generation index
        offset: Width to delay the start of the index
    """
    mapping = reverse(mapping)
    index = [' ' * offset for _ in range(len(str(n)))]
    for number in map(str, (mapping[i] if i in mapping else ' ' for i in range(n))):
        for idx, n in enumerate(number):
            index[idx] += n + '|'
    pad = max(map(len, index))
    return '\n'.join(i.rjust(pad, ' ') for i in reversed(index))


def vertical_index(n, mapping):
    """Create a generator that return the next generation index for that array row

    Args:
        n: Number of vertical indexes to generate
        mappinG: A mapping of top->bottom array index to equivalent generation index"""
    mapping = reverse(mapping)
    width = len(str(n))
    index = map(str, (mapping[i] if i in mapping else ' ' * len(str(i)) for i in range(n)))
    for ind in (i.rjust(width, ' ') + '|' for i in index):
        yield '\n' + ind



