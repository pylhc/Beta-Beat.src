import sys
import os

# Root of Beta-Beat.src
new_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
if new_path not in sys.path:
    sys.path.append(new_path)

new_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../Beta-Beat.src/Python_Classes4MAD"))
if new_path not in sys.path:
    sys.path.append(new_path)


new_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../Beta-Beat.src/"))
if new_path not in sys.path:
    sys.path.append(new_path)

new_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
if new_path not in sys.path:
    sys.path.append(new_path)
    
new_path = os.path.abspath(os.path.join( os.path.dirname(os.path.abspath(__file__)), "..", "..", "RDTs" ))
if new_path not in sys.path:
    sys.path.append(new_path)

