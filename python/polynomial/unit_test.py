#!/usr/bin/python

import unittest
import polynomial as pm

class TestPolynomial(unittest.TestCase):

    def setUp(self):
        self.str = str()
        self.pm1 = pm.Polynomial([1, 2, 3, 4])
        self.pm2 = pm.Polynomial([8.17, 9.34, 0, -17])
        self.pm3 = pm.Polynomial([5, 0, 0, 0, 0, 0, 1])
        pass

    # __str__
    def test_method__str__long(self):
        self.pm1 = pm.Polynomial(range(10))
        self.assertEqual(self.pm1.__str__(), "x^8+2x^7+3x^6+4x^5+5x^4+6x^3+7x^2+8x+9")

    def test_method__str__user_err_str(self):
        self.pm1 = pm.Polynomial("10, 2.33, 3.44,    5.00")
        self.assertEqual(self.pm1.__str__(), "10x^3+2.33x^2+3.44x+5")

    def test_method__str__single_pos_value(self):
        self.pm1 = pm.Polynomial(1000)
        self.assertEqual(self.pm1.__str__(), "1000")

    def test_method__str__single_neg_value(self):
        self.pm1 = pm.Polynomial([-100])
        self.assertEqual(self.pm1.__str__(), "-100")

    def test_method__str__tuple(self):
        self.pm1 = pm.Polynomial((1, 2, 3, 4))
        self.assertEqual(self.pm1.__str__(), "x^3+2x^2+3x+4")

    def test_method__str__list(self):
        self.assertEqual(self.pm2.__str__(), "8.17x^3+9.34x^2-17")

    def test_method__str__neg_values(self):
        self.pm1 = pm.Polynomial([-100, -200, -300, -400])
        self.assertEqual(self.pm1.__str__(), "-100x^3-200x^2-300x-400")

    def test_method__str__long_frac_values(self):
        self.pm1 = pm.Polynomial("-1.23456789, 2.555, 3.496, 4.594, 5.50000001")
        self.assertEqual(self.pm1.__str__(), "-1.23x^4+2.56x^3+3.5x^2+4.59x+5.5")

    # __eq__
    def test_method__eq__long_false(self):
        self.pm1 = pm.Polynomial(range(10))
        self.pm2 = pm.Polynomial(range(11))
        self.assertFalse(self.pm1 == self.pm2, "Error! Polynomials aren't equal")

    def test_method__eq__long_true(self):
        self.pm1 = pm.Polynomial(range(100))
        self.pm2 = pm.Polynomial(range(100))
        self.assertTrue(self.pm1 == self.pm2, "Error! Polynomials are equal")

    def test_method__eq__list_tupl_cmp(self):
        self.pm2 = pm.Polynomial((1, 2, 3, 4))
        self.assertTrue(self.pm1 == self.pm2, "Error! Polynomials are equal")

    def test_method__eq__list_str_cmp(self):
        self.pm2 = pm.Polynomial("1, 2, 3, 4")
        self.assertTrue(self.pm1 == self.pm2, "Error! Polynomials are equal")

    def test_method__eq__list_list_cmp(self):
        self.pm1 = pm.Polynomial([8.17, 9.34, 0, -17])
        self.assertTrue(self.pm1 == self.pm2, "Error! Polynomials are equal")

    def test_method__eq__cmp_different_fracs_1(self):
        self.pm1 = pm.Polynomial([8.17, 1.234, 1.239, 0.00001])
        self.pm2 = pm.Polynomial([8.17, 1.23, 1.24, 0])
        self.assertFalse(self.pm1 == self.pm2, "Error! Polynomials aren't equal")

    def test_method__eq__cmp_same_fracs_2(self):
        self.pm1 = pm.Polynomial([8.17, 1.234, 1.239, 0.00001])
        self.pm2 = pm.Polynomial("8.17, 1.234, 1.239, 0.00001")
        self.assertTrue(self.pm1 == self.pm2, "Error! Polynomials aren't equal")

    def test_method__eq__single_val(self):
        self.pm1 = pm.Polynomial([1])
        self.pm2 = pm.Polynomial("1")
        self.assertTrue(self.pm1 == self.pm2, "Error! Polynomials aren't equal")

    # __ne__
    def test_method__ne__long_false(self):
        self.pm1 = pm.Polynomial(range(10))
        self.pm2 = pm.Polynomial(range(11))
        self.assertTrue(self.pm1 != self.pm2, "Error! Polynomials aren't equal")

    def test_method__ne__long_true(self):
        self.pm1 = pm.Polynomial(range(100))
        self.pm2 = pm.Polynomial(range(100))
        self.assertFalse(self.pm1 != self.pm2, "Error! Polynomials are equal")

    def test_method__ne__list_tupl_cmp(self):
        self.pm2 = pm.Polynomial((1, 2, 3, 4))
        self.assertFalse(self.pm1 != self.pm2, "Error! Polynomials are equal")

    def test_method__ne__list_str_cmp(self):
        self.pm2 = pm.Polynomial("1, 2, 3, 4")
        self.assertFalse(self.pm1 != self.pm2, "Error! Polynomials are equal")

    def test_method__ne__list_list_cmp(self):
        self.pm1 = pm.Polynomial([8.17, 9.34, 0, -17])
        self.assertFalse(self.pm1 != self.pm2, "Error! Polynomials are equal")

    def test_method__ne__cmp_different_fracs_1(self):
        self.pm1 = pm.Polynomial([8.17, 1.234, 1.239, 0.00001])
        self.pm2 = pm.Polynomial([8.17, 1.23, 1.24, 0])
        self.assertTrue(self.pm1 != self.pm2, "Error! Polynomials aren't equal")

    def test_method__ne__cmp_same_fracs_2(self):
        self.pm1 = pm.Polynomial([8.17, 1.234, 1.239, 0.00001])
        self.pm2 = pm.Polynomial("8.17, 1.234, 1.239, 0.00001")
        self.assertFalse(self.pm1 != self.pm2, "Error! Polynomials aren't equal")

    def test_method__ne__single_val(self):
        self.pm1 = pm.Polynomial([1])
        self.pm2 = pm.Polynomial("1")
        self.assertFalse(self.pm1 != self.pm2, "Error! Polynomials aren't equal")

    #__iadd__
    def test_method__iadd__single_val(self):
        self.pm1 = pm.Polynomial([1])
        self.pm2 = pm.Polynomial("1")
        self.pm1 += self.pm2
        self.assertEqual(self.pm1.__str__(), "2")

    def test_method__iadd__different_coeffs_num_1(self):
        self.pm1 = pm.Polynomial([1, 2, 3, 4])
        self.pm2 = pm.Polynomial("1, 2, 3, 4, 5, 6")
        self.pm1 += self.pm2
        self.assertEqual(self.pm1.__str__(), "x^5+2x^4+4x^3+6x^2+8x+10")

    def test_method__iadd__different_coeffs_num_2(self):
        self.pm1 = pm.Polynomial([1, 2, 3, 4, 5, 6])
        self.pm2 = pm.Polynomial("1, 2, 3, 4")
        self.pm1 += self.pm2
        self.assertEqual(self.pm1.__str__(), "x^5+2x^4+4x^3+6x^2+8x+10")

    def test_method__iadd__long_coeffs_num(self):
        self.pm1 = pm.Polynomial("1.00001, 2.999999, 3.44, 2.01, 0.00001, 1")
        self.pm2 = pm.Polynomial("9.00001, 3.0, 1.56, 10.01, 0.00001, 1")
        self.pm1 += self.pm2
        self.assertEqual(self.pm1.__str__(), "10x^5+6x^4+5x^3+12.02x^2+2")

    def test_method__iadd__pm_and_float(self):
        self.pm1 = pm.Polynomial([1, 2, 3, 4, 5, 6])
        self.pm1 += 1
        self.assertEqual(self.pm1.__str__(), "x^5+2x^4+3x^3+4x^2+5x+7")

    #__add__
    def test_method__add__single_val(self):
        self.pm1 = pm.Polynomial([1])
        self.pm2 = pm.Polynomial("1")
        self.pm3 = self.pm1 + self.pm2
        self.assertEqual(self.pm3.__str__(), "2")

    def test_method__add__different_coeffs_num_1(self):
        self.pm1 = pm.Polynomial([1, 2, 3, 4])
        self.pm2 = pm.Polynomial("1, 2, 3, 4, 5, 6")
        self.pm3 = self.pm1 + self.pm2
        self.assertEqual(self.pm3.__str__(), "x^5+2x^4+4x^3+6x^2+8x+10")

    def test_method__add__different_coeffs_num_2(self):
        self.pm1 = pm.Polynomial([1, 2, 3, 4, 5, 6])
        self.pm2 = pm.Polynomial("1, 2, 3, 4")
        self.pm3 = self.pm1 + self.pm2
        self.assertEqual(self.pm3.__str__(), "x^5+2x^4+4x^3+6x^2+8x+10")

    def test_method__add__long_coeffs(self):
        self.pm1 = pm.Polynomial("1.00001, 2.999999, 3.44, 2.01, 0.00001, 1")
        self.pm2 = pm.Polynomial("9.00001, 3.0, 1.56, 10.01, 0.00001, 1")
        self.pm3 = self.pm1 + self.pm2
        self.assertEqual(self.pm3.__str__(), "10x^5+6x^4+5x^3+12.02x^2+2")

    def test_method__add__pm_and_float(self):
        self.pm1 = pm.Polynomial([1, 2, 3, 4, 5, 6])
        self.pm2 = self.pm1 + 2
        self.assertEqual(self.pm2.__str__(), "x^5+2x^4+3x^3+4x^2+5x+8")

    #__isub__
    def test_method__isub__single_val(self):
        self.pm1 = pm.Polynomial([1])
        self.pm2 = pm.Polynomial("1")
        self.pm1 -= self.pm2
        self.assertEqual(self.pm1.__str__(), "0")

    def test_method__isub__different_coeffs_num_1(self):
        self.pm1 = pm.Polynomial([1, 2, 3, 4])
        self.pm2 = pm.Polynomial("1, 2, 3, 4, 5, 6")
        self.pm1 -= self.pm2
        self.assertEqual(self.pm1.__str__(), "-x^5-2x^4-2x^3-2x^2-2x-2")

    def test_method__isub__different_coeffs_num_2(self):
        self.pm1 = pm.Polynomial([1, 2, 3, 4, 5, 6])
        self.pm2 = pm.Polynomial("1, 2, 3, 4")
        self.pm1 -= self.pm2
        self.assertEqual(self.pm1.__str__(), "x^5+2x^4+2x^3+2x^2+2x+2")

    def test_method__isub__long_coeffs_num(self):
        self.pm1 = pm.Polynomial("1.00001, 2.999999, 3.44, 2.01, 0.00001, 1")
        self.pm2 = pm.Polynomial("9.00001, 3.0, 1.56, 10.01, 0.00001, 1")
        self.pm1 -= self.pm2
        self.assertEqual(self.pm1.__str__(), "-8x^5+1.88x^3-8x^2")

    def test_method__isub__pm_and_float(self):
        self.pm1 = pm.Polynomial([1, 2, 3, 4, 5, 6])
        self.pm1 -= 1
        self.assertEqual(self.pm1.__str__(), "x^5+2x^4+3x^3+4x^2+5x+5")

    #__mul__
    def test_method__sub__single_val(self):
        self.pm1 = pm.Polynomial([1])
        self.pm2 = pm.Polynomial("1")
        self.pm3 = self.pm1 - self.pm2
        self.assertEqual(self.pm3.__str__(), "0")

    def test_method__sub__different_coeffs_num_1(self):
        self.pm1 = pm.Polynomial([1, 2, 3, 4])
        self.pm2 = pm.Polynomial("1, 2, 3, 4, 5, 6")
        self.pm3 = self.pm1 - self.pm2
        self.assertEqual(self.pm3.__str__(), "-x^5-2x^4-2x^3-2x^2-2x-2")

    def test_method__sub__different_coeffs_num_2(self):
        self.pm1 = pm.Polynomial([1, 2, 3, 4, 5, 6])
        self.pm2 = pm.Polynomial("1, 2, 3, 4")
        self.pm3 = self.pm1 - self.pm2
        self.assertEqual(self.pm3.__str__(), "x^5+2x^4+2x^3+2x^2+2x+2")

    def test_method__sub__long_coeffs_num(self):
        self.pm1 = pm.Polynomial("1.00001, 2.999999, 3.44, 2.01, 0.00001, 1")
        self.pm2 = pm.Polynomial("9.00001, 3.0, 1.56, 10.01, 0.00001, 1")
        self.pm3 = self.pm1 - self.pm2
        self.assertEqual(self.pm3.__str__(), "-8x^5+1.88x^3-8x^2")

    def test_method__sub__pm_and_float(self):
        self.pm1 = pm.Polynomial([1, 2, 3, 4, 5, 6])
        self.pm2 = self.pm1 - 2
        self.assertEqual(self.pm2.__str__(), "x^5+2x^4+3x^3+4x^2+5x+4")

    #__imul__
    def test_method__imul__single_val(self):
        self.pm1 = pm.Polynomial([1])
        self.pm1 *= pm.Polynomial("1")
        self.assertEqual(self.pm1.__str__(), "1")

    def test_method__imul__different_coeffs_num_1(self):
        self.pm1 = pm.Polynomial((-1, 0, 1))
        self.pm1 *= (3, 2, 0, -1)
        self.assertEqual(self.pm1.__str__(), "-3x^5-2x^4+3x^3+3x^2-1")

    def test_method__imul__different_coeffs_num_2(self):
        self.pm1 = pm.Polynomial((3, 2, 0, -1))
        self.pm1 *= (-1, 0, 1)
        self.assertEqual(self.pm1.__str__(), "-3x^5-2x^4+3x^3+3x^2-1")

    def test_method__imul__long_coeffs_num(self):
        self.pm1 = pm.Polynomial("1.00001, 2.999999, 3.44, 2.01, 0.00001, 1")
        self.pm2 = pm.Polynomial("9.00001, 3.0, 1.56, 10.01, 0.00001, 1")
        self.pm1 *= self.pm2
        self.assertEqual(self.pm1.__str__(), "9x^10+30x^9+41.52x^8+43.1x^7+41.43x^6+47.57x^5+26.12x^4+5x^3+12.02x^2+1")

    def test_method__imul__pm_and_float(self):
        self.pm1 = pm.Polynomial([1, 2, 3, 4, 5, 6])
        self.pm1 *= 2
        self.assertEqual(self.pm1.__str__(), "2x^5+4x^4+6x^3+8x^2+10x+12")

    #__mul__
    def test_method__mul__single_val(self):
        self.pm1 = pm.Polynomial([1])
        self.pm2 = self.pm1 * pm.Polynomial("1")
        self.assertEqual(self.pm2.__str__(), "1")

    def test_method__mul__different_coeffs_num_1(self):
        self.pm1 = pm.Polynomial((-1, 0, 1))
        self.pm2 = self.pm1 * (3, 2, 0, -1)
        self.assertEqual(self.pm2.__str__(), "-3x^5-2x^4+3x^3+3x^2-1")

    def test_method__mul__different_coeffs_num_2(self):
        self.pm1 = pm.Polynomial((3, 2, 0, -1))
        self.pm2 = self.pm1 * (-1, 0, 1)
        self.assertEqual(self.pm2.__str__(), "-3x^5-2x^4+3x^3+3x^2-1")

    def test_method__mul__long_coeffs_num(self):
        self.pm1 = pm.Polynomial("1.00001, 2.999999, 3.44, 2.01, 0.00001, 1")
        self.pm2 = pm.Polynomial("9.00001, 3.0, 1.56, 10.01, 0.00001, 1")
        self.pm3 = self.pm1 * self.pm2
        self.assertEqual(self.pm3.__str__(), "9x^10+30x^9+41.52x^8+43.1x^7+41.43x^6+47.57x^5+26.12x^4+5x^3+12.02x^2+1")

    def test_method__mul__pm_and_float(self):
        self.pm1 = pm.Polynomial([1, 2, 3, 4, 5, 6])
        self.pm2 = self.pm1 * 2
        self.assertEqual(self.pm2.__str__(), "2x^5+4x^4+6x^3+8x^2+10x+12")

    def tearDown(self):
        del self.str
        del self.pm1
        del self.pm2
        del self.pm3
        pass

if __name__ == '__main__':
    unittest.main()
