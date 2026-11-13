"""JSON for large test-bearing .esm documents: indented like json.dumps(indent=2), except that
a list of numbers (at any depth of nesting) and a test assertion are each written on one line.
Same data, a fraction of the lines; used where a document carries whole-column test states
(couplings/mosaic_column.esm).

    from esm_compact_json import dumps
    open(path, 'w').write(dumps(doc) + '\n')
"""
import json


def _numeric(x):
    if isinstance(x, list):
        return all(_numeric(y) for y in x)
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def _inline(x):
    return json.dumps(x, separators=(', ', ': '))


def dumps(x, ind=0):
    pad, pad1 = '  ' * ind, '  ' * (ind + 1)
    if isinstance(x, dict):
        if not x:
            return '{}'
        if 'variable' in x and 'expected' in x:          # a test assertion
            return _inline(x)
        items = [f'{pad1}{json.dumps(k)}: {dumps(v, ind + 1)}' for k, v in x.items()]
        return '{\n' + ',\n'.join(items) + '\n' + pad + '}'
    if isinstance(x, list):
        if not x:
            return '[]'
        if _numeric(x):
            return _inline(x)
        return '[\n' + ',\n'.join(pad1 + dumps(v, ind + 1) for v in x) + '\n' + pad + ']'
    return json.dumps(x)
