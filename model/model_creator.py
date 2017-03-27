from __future__ import print_function
import os
import sys
import logging
import creators_list

CURRENT_DIR = os.path.dirname(__file__)
sys.path.append(os.path.abspath(os.path.join(CURRENT_DIR, "..")))
from madx import madx_wrapper

LOGGER = logging.getLogger(__name__)


def _i_am_main():
    accel = sys.argv[1].lower()
    try:
        creator = creators_list.CREATORS[accel]
    except KeyError:
        raise ModelCreationError(
            "First argument should be one of: " +
            str(creators_list.CREATORS.keys())
        )
    creator.start_from_terminal()


class ModelCreator(object):

    @classmethod
    def create_model(creator, instance, output_path):
        instance.verify_object()
        madx_script = creator.get_madx_script(
            instance,
            output_path
        )
        creator.run_madx(madx_script)

    @staticmethod
    def run_madx(madx_script):
        madx_wrapper.resolve_and_run_string(madx_script)


class ModelCreationError(Exception):
    pass


if __name__ == "__main__":
    _i_am_main()
