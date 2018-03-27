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
        elif (type(_coeffs) is list) or (type(_coeffs) is tuple):
            self.parse_list_coeffs(_coeffs)
        else:
            raise ValueError("Error! Unknown type of coefficients - " + type(coeffs))
        print(self.coeffs)

    def polynomial_to_string(self, var_string = 'x', fraction = 2):
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
                else:
                    int_str = str(str(coeff).split('.')[0])
                    frac_str = str(int(str(coeff).split('.')[1]) % 10 ** fraction)
                    str_coeff = int_str + '.' + frac_str

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
        if len(self) != len(other):
            return False
        for i, j in zip(self.coeffs, other.coeffs):
            if i != j:
                return False
        return True

    

if __name__ == "__main__":
    print("Main")
    print(Polynomial([2.2349, 0, -1, 3]) == Polynomial([2.2349, 0, 1, 3]))
