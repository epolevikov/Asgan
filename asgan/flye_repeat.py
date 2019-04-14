'''
import os


def run_flye_repeat(sequences_fasta, out_dir, out_name):
    FLYE_REPEAT_PY = "lib/Flye/bin/flye-repeat-standalone.py"

    cmd = ["python2.7", FLYE_REPEAT_PY]
    cmd.extend(["--input-assembly", sequences_fasta])
    cmd.extend(["--out-dir", out_dir])
    cmd = " ".join(cmd)

    os.system(cmd)
    os.system("mv {0}/assembly_graph.gfa {0}/{1}".format(out_dir, out_name))

    return "{}/{}".format(out_dir, out_name)
'''
