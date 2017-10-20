import os
import pytest


def test_all():
    pytest.main([os.path.join(
        os.path.dirname(__file__),
        "test",
    )])


if __name__ == "__main__":
    test_all()
