import sys
import os

CURRENT_DIR = os.path.dirname(__file__)
sys.path.append(os.path.abspath(os.path.join(CURRENT_DIR, "..")))

from model_creators import model_creator
from model_creators.lhc_model_creator import (LhcModelCreator,
                                              LhcBestKnowledgeCreator)


CREATORS = {
    "lhc": LhcModelCreator,
    "lhc_best_knowledge": LhcBestKnowledgeCreator,
}


def _i_am_main():
    accel = sys.argv[1].lower()
    try:
        creator = CREATORS[accel]
    except KeyError:
        raise model_creator.ModelCreationError(
            "First argument should be one of: " +
            str(CREATORS.keys())
        )
    creator.start_from_terminal()


if __name__ == "__main__":
    _i_am_main()
