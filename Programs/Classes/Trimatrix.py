class Trimatrix:

    @staticmethod
    def get_size(dimension):
        return int(dimension*(dimension-1)/2)

    @staticmethod
    def get_index(i, j):  # conditions: i!=j and (i, j) < dimension
        return int((i*(i-1)/2) + j) if i > j else int((j*(j-1)/2) + i)
