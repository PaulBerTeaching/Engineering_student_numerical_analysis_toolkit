"""
Created on Fri Apr 24 00:40:30 2020

@author: paulberraute
"""

import pytest
import sys

sys.path.append("..")

from dichotomie import dichotomie

def test_char():
    assert dichotomie(0, 5, 100) == 1.3652300134140969


if __name__ == '__main__':
    try:
        test_char()
    except:
        print("Votre fonction de dichotomie ne retourne pas exactement la valeur attendu. Revérifiez s'il vous plaît.")