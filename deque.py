from typing import Sequence
from collections import deque

Deque = Sequence

def rotated_deque(_deque: Deque, n: int) -> Deque:
    d = deque(_deque)
    d.rotate(n)
    return d

def reversed_deque(_deque: Deque, n: int) -> Deque:
    d = deque(_deque)
    d.reverse()
    return d
