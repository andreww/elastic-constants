# Test the setup code by comparing generated files
# with files previously generated for jadeite. Do 
# comparison on line by line and word by word basis
# allowing for numerical errors in numerical values

import os
import shutil
import subprocess
import glob

import numpy as np

def compare_files(file1, file2):
    """
    Line by line, word by word diff with roundoff

    Reads two files and does nothing if they match,
    raises an exeption if they differ
    """
    with open(file1) as f1:
        with open(file2) as f2:
            for line1, line2 in zip(f1, f2):
                words1 = line1.split()
                words2 = line2.split()
                assert len(words1) == len(words2), "line length:\n" + line1 + "\n" + line2
                for word1, word2 in zip(words1, words2):
                    try:
                        num1 = float(word1)
                        num2 = float(word2)
                        np.testing.assert_allclose(num1, num2)
                    except:
                        assert word1 == word2, "words:\n" + line1 + "\n" + line2


def test_jadeite_setup(tmpdir):
    """
    Run generate_strain.py with the jadeite input

    In a tempdir run generate strain and check the files
    that are made match (to high precission) reference 
    files.
    """
    tmpdir = str(tmpdir) # We may run in python 2 with old os
    input_path = os.path.dirname(os.path.realpath(__file__)) + "/generate_strain_files/Input"
    reference_path = os.path.dirname(os.path.realpath(__file__)) + "/generate_strain_files/Output"
    shutil.copytree(input_path, tmpdir + "/input")

    # Run analysis and check results
    olddir = os.getcwd()
    try:
        os.chdir(tmpdir + "/input")
        return_code = subprocess.call([olddir + "/generate_strain.py", "-s 0.1", "jadeite_4_30_GPa"])

        files = glob.glob("*")
        for file in files:
            compare_files(file, reference_path + "/" + file)

    finally:
        os.chdir(olddir)


if __name__ == "__main__":
    import sys
    compare_files(sys.argv[1], sys.argv[2])
