from __future__ import print_function
import logging
from madx import madx_wrapper

LOGGER = logging.getLogger(__name__)


class ModelCreator(object):

    @classmethod
    def create_model(creator, instance, output_path):
        instance.verify_object()
        madx_script = creator.get_madx_script(
            instance,
            output_path
        )
        creator.prepare_run(instance, output_path)
        creator.run_madx(madx_script)

    @staticmethod
    def run_madx(madx_script):
        madx_wrapper.resolve_and_run_string(madx_script)


class ModelCreationError(Exception):
    pass
