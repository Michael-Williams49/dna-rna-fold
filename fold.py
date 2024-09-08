# Define the input sequence
sequence = "GCAGCUGCCAUCUUAGGGGCGCCUGGCGCUACGGGUUUCUCGUUGGAGGCGGCCUUCGUGGCAGCUGUAGACGCCGGGAAAAGGCAUAAAGUCCGUUGGCCGAC"

# Define the scoring matrix for base pair matching
match = {
    "A": {"A": -5, "G": -5, "C": -5, "T": 2, "U": 2},
    "G": {"A": -5, "G": -5, "C": 3, "T": -5, "U": 1},
    "C": {"A": -5, "G": 3, "C": -5, "T": -5, "U": -5},
    "T": {"A": 2, "G": -5, "C": -5, "T": -5, "U": -5},
    "U": {"A": 2, "G": 1, "C": -5, "T": -5, "U": -5}
}

# Define the penalty for insertions and deletions (indels)
indel = {
    False: -5,  # Penalty for non-consecutive indels
    True: -1    # Penalty for extension of indels
}

# Define the minimum score threshold for considering a base pair
threshold = 7

# Define the Scheme class to encapsulate scoring parameters
class Scheme:
    def __init__(self, match, indel, threshold):
        self.match = match
        self.indel = indel
        self.threshold = threshold

# Define the Table class for storing and manipulating the dynamic programming matrix
class Table:
    def __init__(self, sequence, default):
        self.sequence = sequence
        self.inverse = sequence[::-1]  # Reverse sequence for complementary base pairing
        self.default = default
        self.length = len(sequence)
        # Initialize the 2D table with default values
        self.data = {i: {j: default for j in range(self.length - i + 1)} for i in range(self.length + 1)}

    def __getitem__(self, index):
        return self.data[index]
    
    def __repr__(self):
        # Create a string representation of the table for debugging
        representation = ""
        representation += "\t-"
        for i in range(self.length):
            representation += "\t" + self.inverse[i]
        representation += "\n"
        for i in range(self.length + 1):
            if i == 0:
                representation += "-"
            else:
                representation += self.sequence[i - 1]
            for j in range(self.length - i + 1):
                representation += "\t" + str(self.data[i][j])
            representation += "\n"
        return representation
    
    def max(self):
        # Find the maximum value in the table and its position
        maximum = {"value": self.default, "x": 0, "y": 0}
        for i in range(1, self.length + 1):
            for j in range(1, self.length - i + 1):
                if self.data[i][j] > maximum["value"]:
                    maximum["value"] = self.data[i][j]
                    maximum["x"] = i
                    maximum["y"] = j
        return maximum

# Define the Fold class for performing the folding algorithm
class Fold:
    def __init__(self, sequence, scheme):
        self.sequence = sequence
        self.inverse = sequence[::-1]
        self.length = len(sequence)
        self.scheme = scheme
        self.score = Table(sequence, 0)  # Table for storing scores
        self.transition = Table(sequence, 0)  # Table for storing transition information
        self.pair = []  # List to store paired bases
        self.mask = {"x": {0}, "y": {0}}  # Set of masked positions
        self.control()  # Start the folding process

    def __repr__(self):
        # Create a string representation of the folding result
        path = [0 for i in range(self.length + 1)]
        for point in self.pair:
            path[point[0]] = point[1]
            path[point[1]] = point[0]
        representation = ""
        for i in range(self.length):
            representation += str(i + 1) + "\t"
        representation += "\n"
        for i in range(self.length):
            representation += self.sequence[i] + "\t"
        representation += "\n"
        for i in range(self.length):
            representation += str(path[i + 1]) + "\t"
        representation += "\n"
        return representation
    
    def fasta(self):
        # Generate a FASTA-like representation of the folding result
        representation = "." * self.length
        for i, j in self.pair:
            representation = representation[:i-1] + "(" + representation[i:]
            representation = representation[:j-1] + ")" + representation[j:]
        return ">sequence\n" + self.sequence + "\n" + representation 
        
    def fill(self):
        # Fill the score and transition tables using dynamic programming
        for i in range(1, self.length + 1):
            for j in range(1, self.length - i + 1):
                if i in self.mask["x"] or j in self.mask["y"]:
                    continue
                values = [
                    self.score[i - 1][j - 1] + self.scheme.match[self.sequence[i - 1]][self.inverse[j - 1]],
                    self.score[i - 1][j] + self.scheme.indel[self.transition[i - 1][j] != 1],
                    self.score[i][j - 1] + self.scheme.indel[self.transition[i][j - 1] != 1]
                ]
                optimum = max(values)
                if optimum <= 0:
                    self.score[i][j] = 0
                    self.transition[i][j] = 0
                else:
                    self.transition[i][j] = [index + 1 for index, value in enumerate(values) if value == optimum][0] # match = 1, down = 2, right = 3
                    self.score[i][j] = optimum

    def trace(self, x, y):
        # Trace back the optimal path from a given position
        path = []
        while not self.score[x][y] == 0:
            transit = self.transition[x][y]
            if transit == 1:
                path.append((x, y))
                x = x - 1
                y = y - 1
            elif transit == 2:
                x = x - 1
            elif transit == 3:
                y = y - 1
        return path
    
    def eliminate(self, x, y):
        # Mask positions that have been paired or should be excluded
        for i in range(1, self.length - y + 1):
            self.score[i][y] = 0
            self.transition[i][y] = 0
        for j in range(1, self.length - x + 1):
            self.score[x][j] = 0
            self.transition[x][j] = 0
        for i in range(1, x):
            self.score[i][self.length - x + 1] = 0
            self.transition[i][self.length - x + 1] = 0
        for j in range(1, y):
            self.score[self.length - y + 1][j] = 0
            self.transition[self.length - y + 1][j] = 0
        self.mask["x"].update([x, self.length - y + 1])
        self.mask["y"].update([y, self.length - x + 1])

    def control(self):
        # Main control loop for the folding algorithm
        self.fill()
        record, x, y = self.score.max().values()
        while record > self.scheme.threshold:
            path = self.trace(x, y)
            self.eliminate(x, y)
            for point in path:
                self.eliminate(*point)
                self.pair.append((point[0], self.length - point[1] + 1))
            self.fill()
            record, x, y = self.score.max().values()
            
# Create a Scheme object with the defined parameters
scheme = Scheme(match, indel, threshold)

# Create a Fold object and perform the folding
fold = Fold(sequence, scheme)

# Print the result in FASTA-like format
print(fold.fasta())