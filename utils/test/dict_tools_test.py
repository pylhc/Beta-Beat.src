import sys
import os
from os.path import abspath, join, dirname
new_path = abspath(join(dirname(abspath(__file__)), os.pardir, os.pardir))
if new_path not in sys.path:
    sys.path.append(new_path)

from utils.dict_tools import DotDict
from utils.dict_tools import DictParser
from utils.dict_tools import Argument
from utils.dict_tools import print_dict_tree


def simple_dictparse_test():
    parser = DictParser({
        'test': Argument(
            name='test',
            default='test',
            required=False,
        ),
        'test1': Argument(
            name='test1'
        ),
        'sub': {
            'test': Argument(
                'test',
                type=int,
                default=2
            )
        },
        'sub2': {
            'subsub': {
                'buub': Argument(
                    'buub',
                    type=int,
                    default=2
                )
            }

        }
    }
    )

    # parser.tree()
    opt = {
        'sub': {
            'test': 2
        },
        'sub2': {
            'subsub': {
                'buub': 2
            }
        }
    }

    opt = parser.parse_options(opt)
    print_dict_tree(opt, "Options")


if __name__ == "__main__":
    simple_dictparse_test()
