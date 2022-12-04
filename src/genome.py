from __future__ import annotations

from abc import (
    # A tag that says that we can't use this class except by specialising it
    ABC,
    # A tag that says that this method must be implemented by a child class
    abstractmethod
)

class Genome(ABC):
    """Representation of a circular genome."""

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

class ListGenome(Genome):
    """
    Representation of a genome.

    Implements the Genome interface using Python's built-in lists
    """
    def __init__(self, n: int):
        """Create a new genome with length n."""
        gen_seq=[] #Making the list that will become the genome list. Its an empty list to make sure the length is n and not n+1.
        for _ in range(n):
            gen_seq+=['-'] #To make the list actually longer i add a string. specifically '-' for the string representation of non-TE_regions
        self.genome = gen_seq
        self.te_dict = {}
        self.te_count = 0

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
        self.te_count += 1 #If we are inserting a TE the count/ID should change
        self.te_dict[self.te_count] = [pos, length] #Insertion of the te into our dict with position as key and length as value
        te_positions = list(self.te_dict.keys()) #Making a list of keys for identifying if the te to be inserted overlap with 
        # an exsisting/active te.
        for position in te_positions: #itterate over each element in the te dictionary
            start = self.te_dict[position][0] # defining a start position
            end = self.te_dict[position][1] + start # defining the end position to know the range
        #    print(position,self.te_dict, 'for1',start,pos,end)
            if start < pos <= end: # cheking overlap.
        #        print(position, self.te_dict,self.te_dict,'if1')
        #                    print(position, self.te_dict,'while')
                self.disable_te(position) #Disabling the alredy exsiting te (since the sequence get changed)
        #        print(position, self.te_dict,'for2 after disable')
            if start > pos: #The positions are start positions so incase start is bigger than pos we need to update the possitions:
                self.te_dict[position][0] = start + length #adding the length of the inserted TE to the start position, and updating in the dict.
        self.genome[pos:pos] = length*['A'] #Adding the TE to our genome map.
        #print(position, te_positions, 'end')
        #self.disable_te(position) #Disabling the alredy exsiting te (since the sequence get changed)
        #print(position, self.te_dict,pos, length, 'end')
        #print(''.join(self.genome), pos, length, self.te_count, te_positions, position)
        return self.te_count

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
        if te not in list(self.te_dict.keys()): #Making sure the te input exsists as an active te in the genome.
            return None
        else: #using else statement to avoide copying incase its not in the genome
            if offset < 0:
                pos_after_offset = 50 + self.te_dict[te][0] + offset #adding (so either negativ or positive should work) the offset to the position pre offset
            #    print(self.te_dict, self.te_dict[te], self.te_dict[te][0],offset,pos_after_offset,self.te_dict[te][1], "HERE!")
                return self.insert_te(pos_after_offset, self.te_dict[te][1])    
            pos_after_offset = self.te_dict[te][0] + offset #adding (so either negativ or positive should work) the offset to the position pre offset
            #print(self.te_dict, self.te_dict[te], self.te_dict[te][0],offset,pos_after_offset,self.te_dict[te][1], "HERE!")
            return self.insert_te(pos_after_offset, self.te_dict[te][1]) #we want an int output which i assume to be the new id so i just throw it to insert_tee.


    def disable_te(self, te: int) -> None:
        """
        Disable a TE.

        If te is an active TE, then make it inactive. Inactive
        TEs are already inactive, so there is no need to do anything
        for those.
        """
        #print(''.join(self.genome), 'DISABLE!', self.te_dict)
        for nucleotide in range(self.te_dict[te][0], self.te_dict[te][0]+self.te_dict[te][1]): #iterating through all nucleotides in the now inactive TE
                    if self.genome[nucleotide] == 'A':
        #                print(position, self.te_dict,'if2')
                        self.genome[nucleotide] = 'x'
                    if nucleotide +1 == self.te_dict[te][0]+self.te_dict[te][1]:
        #                print(position, self.te_dict,'if3')
                        i = 1
                        while self.genome[nucleotide+i] == 'A':
                            self.genome[nucleotide+i] = 'x'
                            i+=1
        self.te_dict.pop(te)
        #print(self.te_dict)
        return None

    def active_tes(self) -> list[int]:
        """Get the active TE IDs."""
        return list(self.te_dict.keys())

    def __len__(self) -> int:
        """Current length of the genome."""
        return len(Genome)

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
        print(''.join(self.genome))
        return ''.join(self.genome)

#-------------------------------------------Linked list genome representation------------------------------------------------------------#
from typing import(
    Generic, TypeVar
)

T = TypeVar('T')

class Link(Generic[T]):
    """Doubly linked link."""

    val: T
    prev: Link[T]
    next: Link[T]

    def __init__(self, val: T, p: Link[T], n: Link[T]):
        """Create a new link and link up prev and next."""
        self.val = val
        self.prev = p
        self.next = n


def insert_after(link: Link[T], val: T) -> None:
    """Add a new link containing val after link."""
    new_link = Link(val, link, link.next)
    new_link.prev.next = new_link
    new_link.next.prev = new_link


def remove_link(link: Link[T]) -> None:
    """Remove link from the list."""
    link.prev.next = link.next
    link.next.prev = link.prev

class LinkedListGenome(Genome):
    """
    Representation of a genome.

    Implements the Genome interface using linked lists.
    """

    head = Link[T]

    def __init__(self, n: int):
        """Create a new genome with length n."""
        self.head = Link(None, None, None)
        self.head.next = self.head
        self.head.prev = self.head
        self.length = n
        self.te_dict = {}
        self.te_count = -1

            for _ in range(self.length):
                insert_after(self.head.prev, '-')

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
        active_tes = list(self.te_identities.items())
        for id, te_range in active_tes:
            start, end = te_range
        self.te_count += 1
        self.te_dict[self.te_count] = [pos, length]
        te_positions = list(self.te_dict.keys())
        for position in te_positions:
            start = self.te_dict[position][0]
            end = self.te_dict[position][1] + start
            if pos > (start - 1) and pos < (end + 1):
                self.disable_te(position)
            elif start > pos:
                self.te_dict[position] = (start + length, end + length)
        self.te_dict[self.te_count] = (pos, pos + length)
        link1 = self.head.next
        for _ in range(pos - 1):
            link1 = link1.next
        for _ in range(length):
            insert_after(link1, 'A')
            link1 = link1.next
        self.length += length

        return self.te_dict[self.te_count][0]

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
        if te not in self.te_dict:
            return None
        else:
            start = self.te_dict[te][0]
            end = self.te_dict[te][1] + start
            if start + offset < 0:
                return self.insert_te(start + offset + self.length, end - start)
            else:
                return self.insert_te(start + offset, end - start)

    def disable_te(self, te: int) -> None:
        """
        Disable a TE.

        If te is an active TE, then make it inactive. Inactive
        TEs are already inactive, so there is no need to do anything
        for those.
        """
        for _ in range(self.te_dict[te][0]):
            self.head.next = self.head.next.next
        for _ in range((self.te_dict[te][0] + self.te_dict[te][1]) - self.te_dict[te][0]):
            self.head.next.val = 'x'
            self.head.next = self.head.next.next

    def active_tes(self) -> list[int]:
        """Get the active TE IDs."""
        return list(self.te_dict.keys())

    def __len__(self) -> int:
        """Current length of the genome."""
        return self.length

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
        genome_list = []
        for _ in range(self.length):
            genome_list.append(self.head.next.val)
            self.head.next = self.head.next.next
        return "".join(genome_list)
