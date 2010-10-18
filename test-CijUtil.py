import CijUtil
import numpy as np
import unittest

class TestInvertCijFunctions(unittest.TestCase):

    def setUp(self):
        self.inmatrix = np.matrix([[0.700, 0.200],[0.400, 0.600]])
        self.inerrors = np.matrix([[0.007, 0.002],[0.004, 0.006]])
        self.true_inv = np.matrix([[1.765, -0.588],[-1.177, 2.059]])
        self.true_err = np.sqrt(np.matrix([[5.269E-4, 1.603E-4],[6.413E-4, 7.172E-4]]))
        self.true_cov = np.array([[[[5.269E-4,-2.245E-4],[-4.490E-4, 2.514E-4]],
                                   [[-2.245E-4,1.603E-4],[2.514E-4,-2.619E-4]]],
                                  [[[-4.490E-4, 2.514E-4],[6.413E-4, -5.238E-4]],
                                   [[2.514E-4, -2.619E-4],[-5.238E-4,7.172E-4]]]])
        (self.calc_inv, self.calc_err, self.calc_cov) = CijUtil.invertCij(self.inmatrix, self.inerrors)

    def test_inverse(self):
        for i in range(2):
            for j in range(2):
                self.assertAlmostEqual(self.calc_inv[i,j], self.true_inv[i,j], 2)

    def test_inverseErrors(self):
        for i in range(2):
            for j in range(2):
                self.assertAlmostEqual(self.calc_err[i,j], self.true_err[i,j], 4)

    def test_inverseCovar(self):
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    for l in range(2):
                        self.assertAlmostEqual(self.calc_cov[i,j,k,l], self.true_cov[i,j,k,l], 7)

if __name__ == '__main__':
    unittest.main()

