[![Open in Visual Studio Code](https://classroom.github.com/assets/open-in-vscode-c66648af7eb3fe8bc4f294546bfd86ef473780cde1dea487d3c4ff354943c9ae.svg)](https://classroom.github.com/online_ide?assignment_repo_id=9524228&assignment_repo_type=AssignmentRepo)
# Project 5: Simulating transposable elements

In the last project, we imagine that someone has hired us to help out with simulating a genome containing [transposable elements]. (I know people who has such strange interests, so it is not beyond the realm of possibilities).

We won’t do anything complicated, this is just an exercise after all, but we will want to simulate TEs as stretches of DNA that can copy themselves elsewhere in the genome.

Our employer already has most of the simulator up and running. She has a program that randomly picks operations to do—insert a TE ab initio, copy a TE, or disable one with a mutation—but she needs us to program a representation of a genome to track where the TEs are.

There are multiple ways to do this, but you should implement at least two: one based Python lists, where each nucleotide is represented by one entry in a list, and one based on linked lists, where each nucleotide is represented by a link. If you feel ambitious, you can try others (for example keeping track of ranges of a genome with the same annotation so you don’t need to explicitly represent each nucleotide).

## Genome interface

A genome should be represented as a class that implements the following methods:

```python
class Genome(ABC):
    """Representation of a circular enome."""

    def __init__(self, n: int):
        """Create a genome of size n."""
        ...  # not implemented yet

    @abstractmethod
    def insert_te(self, pos: int, length: int) -> int:
        """
        Insert a new transposable element.

        Insert a new transposable element at position pos and len
        nucleotide forward.

        If the TE collides with an existing TE, i.e. genome[pos]
        already contains TEs, then that TE should be disabled and
        removed from the set of active TEs.

        Returns a new ID for the transposable element.
        """
        ...  # not implemented yet

    @abstractmethod
    def copy_te(self, te: int, offset: int) -> int | None:
        """
        Copy a transposable element.

        Copy the transposable element te to an offset from its current
        location.

        The offset can be positive or negative; if positive the te is copied
        upwards and if negative it is copied downwards. If the offset moves
        the copy left of index 0 or right of the largest index, it should
        wrap around, since the genome is circular.

        If te is not active, return None (and do not copy it).
        """
        ...  # not implemented yet

    @abstractmethod
    def disable_te(self, te: int) -> None:
        """
        Disable a TE.

        If te is an active TE, then make it inactive. Inactive
        TEs are already inactive, so there is no need to do anything
        for those.
        """
        ...  # not implemented yet

    @abstractmethod
    def active_tes(self) -> list[int]:
        """Get the active TE IDs."""
        ...  # not implemented yet

    @abstractmethod
    def __len__(self) -> int:
        """Get the current length of the genome."""
        ...  # not implemented yet

    @abstractmethod
    def __str__(self) -> str:
        """
        Return a string representation of the genome.

        Create a string that represents the genome. By nature, it will be
        linear, but imagine that the last character is immidiatetly followed
        by the first.

        The genome should start at position 0. Locations with no TE should be
        represented with the character '-', active TEs with 'A', and disabled
        TEs with 'x'.
        """
        ...  # not implemented yet

```

The `ABC` and `@abstractmethod` just means that this class is not something you can use by itself, but that another class must implement the details. In `src/genome.py` you will find templates for a Python list tand a linked list implementation (without the actual implementation, because you have to implement them).

You are free to implement the genome classes however you want, and using whateer auxilary data structures you desire, as long as one uses a Python list with an element for each nucleotide and the other a linked list with a link for each nucleotide. If you want to implement a third (or fourth or fifth...) version, you are very welcome to do so as well.

## Complexity

When you have implemented the two (or more) classes, describe the complexity of each operation as a function of the genome size (at the time of the operation), and the size of the TE involved (and when copying, the offset you are copying). Put the description here:

Notation:
n = genome size
m = TE size

Listegenome operations:
    def __init__(self, n: int):
        gen_seq=[]
        for _ in range(n):
            gen_seq+=['-']
        self.genome = gen_seq
        self.te_dict = {}
        self.te_count = 0
    Complexity: O(n), for vi laver en liste ved at tilføje ['-'] til listen gen_seq n gange. De resterende operationer er O(k) så
    de påvirker ikke complexiteten.

    def insert_te(self, pos: int, length: int) -> int:
        self.te_count += 1
        self.te_dict[self.te_count] = [pos, length]
        te_positions = list(self.te_dict.keys())
        
        for position in te_positions:
            start = self.te_dict[position][0]
            end = self.te_dict[position][1] + start
            if start < pos <= end:
                self.disable_te(position)
            if start > pos:
                self.te_dict[position][0] = start + length

        self.genome[pos:pos] = length*['A']
        return self.te_count
    Complexity: O(n+m), jeg ignorer lige at vi kalder disable, så den primære complexitet er baseret ud fra selve indsætelsen af TE i genomet. Resten har complexiteten O(k) så igen ignoreres de. O(n+m) kommer så fra "'self.genome[pos:pos] = length*['A']'" i det vi putter hele TE'en ind som tager O(m) at indsætte også O(n) for at rygge alle de følgende elementer.

    def copy_te(self, te: int, offset: int) -> int | None:
        if te not in list(self.te_dict.keys()):
            return None
        else:
            if offset < 0:
                pos_after_offset = len(self.genome) + self.te_dict[te][0] + offset
                return self.insert_te(pos_after_offset, self.te_dict[te][1])    
            pos_after_offset = (self.te_dict[te][0] + offset) % len(self.genome)
            return self.insert_te(pos_after_offset, self.te_dict[te][1])
    Complexity: O(k), igen ignoreres at vi kalder en anden funktion (insert_te). Alt ande i copy er beregninger som køre i O(k).

    def disable_te(self, te: int) -> None:
        pos = self.te_dict[te][0]
        if self.te_dict[te][0] > len(self.genome):
            pos = self.te_dict[te][0]
            pos = pos % len(self.genome)
        for nucleotide in range(pos, pos+self.te_dict[te][1]): #iterating through all nucleotides in the now inactive TE
            if self.genome[nucleotide] == 'A':
                self.genome[nucleotide] = 'x'
            if nucleotide +1 == pos+self.te_dict[te][1]:
                i = 0
                while self.genome[nucleotide+i] == 'A':
                    self.genome[nucleotide+i] = 'x'
                    i+=1
        self.te_dict.pop(te)
        return None
    Complexity: O(m), Primært kører alt igen i O(k) men når 'A' i genomet ændres til 'x' løber vi over m (længden af TE'en der deaktiveres) så O(m).

    def active_tes(self) -> list[int]:
        return list(self.te_dict.keys())
    Complexity: O(k), fordi 'list()' kører i O(k), det samme gælder 'keys()'

    def __len__(self) -> int:
        return len(self.genome)
    Complexity: O(k) fordi 'len()' kører i O(k)

    def __str__(self) -> str:
        return ''.join(self.genome)
    Complexity: O(n), fordi 'join()' kører i O(n) fordi den skal køre over hele genomet også tilføje dem til den tomme string, det tager O(k).


    

In `src/simulate.py` you will find a program that can run simulations and tell you actual time it takes to simulate with different implementations. You can use it to test your analysis. You can modify the parameters to the simulator if you want to explore how they affect the running time.
