import sys
from pathlib import Path

ROOT_DIR = Path.cwd().resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR))


