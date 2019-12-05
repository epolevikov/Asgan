import os
import sys

asgan_root = os.path.dirname(os.path.realpath(__file__))
src = os.path.join(asgan_root, "src")
sys.path.insert(0, src)

from asgan.main import main

main()
