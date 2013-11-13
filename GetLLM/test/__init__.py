import sys
import os


# Add directory to the python search path - needed to run script from command line
# Otherwise ScriptRunner won't be found
new_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
if new_path not in sys.path:
    sys.path.append(new_path)

new_path = os.path.abspath(os.path.dirname(os.path.abspath(__file__)))
if new_path not in sys.path:
    sys.path.append(new_path)


