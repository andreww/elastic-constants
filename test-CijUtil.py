import CijUtil
import numpy as np
import numpy.testing as npt
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
        npt.assert_array_almost_equal(self.calc_inv, self.true_inv, 3,
           'Calculated inverse of test matrix is wrong')

    def test_inverseErrors(self):
        npt.assert_array_almost_equal(self.calc_err, self.true_err, 5,
           'Calculated propogated std. errors of test matrix are wrong')

    def test_inverseCovar(self):
        npt.assert_array_almost_equal(self.calc_cov, self.true_cov,7,
           'Calculated propogated var-covar matrix of test matrix is wrong')

if __name__ == '__main__':
    unittest.main()

