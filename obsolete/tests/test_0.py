#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""unit test for brixs."""

# standard imports
import numpy as np
from pathlib import Path

import brixs as br

# backpack
import unittest




class TestSum(unittest.TestCase):
    def test_list_int(self):
        """
        Test that it can sum a list of integers
        """
        data = [1, 2, 3]
        result = sum(data)
        self.assertEqual(result, 6)

if __name__ == '__main__':
    unittest.main()
