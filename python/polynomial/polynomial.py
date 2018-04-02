#!/usr/bin/python

class Polynomial:

    def __init__(self, _coeffs = 0):
        self.coeffs = list()
        self.parse_coeffs(_coeffs)

    def parse_string_coeffs(self, _coeffs_str):
        self.coeffs = [float(i) for i in _coeffs_str.replace(" ", "").split(",")]

    def parse_list_coeffs(self, _coeffs_list):
        for it in _coeffs_list:
            self.coeffs.append(float(it))

    def parse_coeffs(self, _coeffs):
        if type(_coeffs) is str:
            self.parse_string_coeffs(_coeffs)
        elif (type(_coeffs) is int) or (type(_coeffs) is float):
            self.coeffs.append(float(_coeffs))
        elif (type(_coeffs) is list) or (type(_coeffs) is tuple) or (type(_coeffs) is range):
            self.parse_list_coeffs(_coeffs)
        elif (isinstance(_coeffs, Polynomial)):
            self.coeffs.extend(_coeffs.coeffs)
        else:
            raise ValueError("Error! Unknown type of coefficients - " + str(type(_coeffs)))

    def polynomial_to_string(self, var_string = 'x', fraction = 2):
        if all(coeff == 0 for coeff in self.coeffs):
            return '0'

        res_str = ''
        first_pow = len(self) - 1

        for i, coeff in enumerate(self.coeffs):
            power = first_pow - i
            coeff = round(coeff, fraction)

            if coeff:
                if coeff < 0:
                    sign = '-'
                elif coeff > 0:
                    sign = ('+' if res_str else '')
                coeff = abs(coeff)

                if (coeff == 1) and (power != 0):
                    str_coeff = ''
                elif (coeff == 0) and (power == 0):
                    str_coeff = '0'
                else:
                    int_str = str(str(coeff).split('.')[0])
                    frac_str = str(str(coeff - (int(coeff) % int(int_str))).split('.')[1])
                    str_coeff = int_str + ('' if frac_str == '0' else ('.' + frac_str))

                if power == 0:
                    str_power = ''
                elif power == 1:
                    str_power = var_string
                else:
                    str_power = var_string + '^' + str(power)

                res_str += sign + str_coeff + str_power
        return res_str

    def __str__(self):
        return self.polynomial_to_string()

    def __len__(self):
        return len(self.coeffs)

    def __eq__(self, other):
        other = Polynomial(other)
        this = Polynomial(self.coeffs)
        this.extend_coeffs(other)

        for i, j in zip(this.coeffs, other.coeffs):
            if i != j:
                return False
        return True

    def __ne__(self, other):
        other = Polynomial(other)
        return False if self.__eq__(other) == True else True

    def extend_coeffs(self, other):
        if len(self) < len(other):
            self.coeffs.reverse()
            self.coeffs.extend([0] * (len(other) - len(self)))
            self.coeffs.reverse()
        else:
            other.coeffs.reverse()
            other.coeffs.extend([0] * (len(self) - len(other)))
            other.coeffs.reverse()

    # Add operation
    def __iadd__(self, other):
        other = Polynomial(other)
        result_len = max(len(self) , len(other))
        self.extend_coeffs(other)
        for it in range(0, result_len):
            self.coeffs[it] += other.coeffs[it]
        return self

    def __add__(self, other):
        other = Polynomial(other)
        res = Polynomial(self.coeffs)
        res += other
        return res

    def __radd__(self, other):
        other = Polynomial(other)
        res = Polynomial(self.coeffs)
        other += res
        return other

    # Sub operation
    def __isub__(self, other):
        other = Polynomial(other)
        result_len = max(len(self) , len(other))
        self.extend_coeffs(other)
        for it in range(0, result_len):
            self.coeffs[it] -= other.coeffs[it]
        return self

    def __sub__(self, other):
        other = Polynomial(other)
        res = Polynomial(self.coeffs)
        res -= other
        return res

    def __rsub__(self, other):
        other = Polynomial(other.coeffs)
        res = Polynomial(self.coeffs)
        other -= res
        return other

    # Mul operatiuon
    def __imul__(self, other):
        other = Polynomial(other)
        res = Polynomial([0]*(len(other)+len(self)-1))
        for i1, val1 in enumerate(self.coeffs):
            for i2, val2 in enumerate(other.coeffs):
                res.coeffs[i1 + i2] += val1 * val2
        self.coeffs = []
        for it in res.coeffs:
            self.coeffs.append(it)
        return self

    def __mul__(self, other):
        other = Polynomial(other)
        res = Polynomial(self.coeffs)
        res *= other
        return res

    def __rmul__(self, other):
        other = Polynomial(other)
        res = Polynomial(self.coeffs)
        other *= res
        return other

if __name__ == "__main__":
    print("Main")
