# Verify that MDS+ connection is working and that we can retrieve data for a given shot number

import HIREXSR_py as pkg


# simple test
def test_mds_connection():
    final_shotno = 1160930043
    conn = pkg.openTree(final_shotno)

    assert conn is not None, "Failed to connect to MDS+ server"
    assert pkg.currentShot(conn) == final_shotno, "Failed to draw data from MDS+ server"
