import __init__
from Utilities.dict_tools import DotDict
from Utilities.dict_tools import DictParser
from Utilities.dict_tools import Argument
from Utilities.dict_tools import print_dict_tree


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
