import sys
import os

# Root of Beta-Beat.src
new_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
if new_path not in sys.path:
    sys.path.append(new_path)