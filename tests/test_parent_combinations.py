from analysis import parent_combinations


def test_parent_combinations():
    assert parent_combinations(['xx', 'xX'], ['xy', 'xY']) == \
           [('xx', 'xy'), ('xx', 'xY'), ('xX', 'xy'), ('xX', 'xY')]
