from silver.expression_profile import ExpressionProfile
import pandas as pd
import unittest


class TestExpressionProfile(unittest.TestCase):
    def setUp(self):
        # Address of a toy expression profile used for testing
        address = 'test/data/dummy_expression.txt'
        # Read the DataFrame containing the expression profile
        self.data = pd.read_csv(address, sep='\t', index_col='ID')
        # Instantiate an ExpressionProfile object
        self.profile = ExpressionProfile(self.data)
        # The list of expected items in the self.profile object
        self.items = items = [('a', [1.12, 1.13, 0.97, 0.91, 1.05, 1.05,
                                     1.08, 0.98, 0.79, 1.01, 1.05, 0.96]),
                              ('b', [2.09, 2.03, 2.2, 1.99, 2.08, 2.0,
                                     2.04, 1.89, 1.93, 1.95, 1.97, 2.18]),
                              ('c', [3.03, 3.12, 2.82, 3.07, 2.96, 2.95,
                                     2.99, 3.14, 3.16, 3.08, 3.03, 2.82]),
                              ('d', [4.05, 4.0, 3.95, 4.13, 3.95, 4.14,
                                     3.99, 3.92, 4.16, 4.23, 3.86, 3.93]),
                              ('e', [4.94, 4.87, 5.06, 5.07, 4.99, 5.15,
                                     5.22, 4.91, 4.91, 5.01, 5.04, 5.09])]
        # List of the column names for the expression profile
        self.col_names = list('ABCDEFGHIJKL')
        # list of the row names for the expression profile
        self.row_names = list('abcde')
        # Number of columns, i.e. samples, in the expression profile
        self.NUM_COLUMNS = 12

    def test_keys(self):
        # Check if keys() return an iterator
        keys = [k for k, _ in self.items]
        self.assertTrue(hasattr(self.profile.keys(), '__iter__') and
                        hasattr(self.profile.keys(), '__next__'))
        # Check if keys() return correct values
        self.assertListEqual(list(self.profile.keys()), keys)

    def test_shape(self):
        self.assertSequenceEqual(self.profile.shape, (5, 12))

    def test__len__(self):
        size = len(self.items)
        self.assertEqual(len(self.data), size)

    def test_items(self):
        # Check if items returns the correct values
        self.assertListEqual(list(self.profile.items()), self.items)

    def test_values(self):
        values = [value for _, value in self.items]
        self.assertListEqual(list(self.profile.values()), values)

    def test__getitem__(self):
        # Check if the method returns TypeError for an invalid index type.
        with self.assertRaises(TypeError, msg='i must be of type integer.'):
            self.profile['a']
        with self.assertRaises(TypeError, msg='i must be of type integer.'):
            self.profile[2, 3]
        with self.assertRaises(TypeError, msg='i must be of type integer.'):
            self.profile[3.12]
        # Check if the method returns IndexError for out of bound indices.
        with self.assertRaises(IndexError, msg='Index out of bound error.'):
            self.profile[-1]
        with self.assertRaises(IndexError, msg='Index out of bound error.'):
            self.profile[5]
        with self.assertRaises(IndexError, msg='Index out of bound error.'):
            self.profile[6]
        # Check if the method returns the right value for valid indices.
        for i in range(len(self.items)):
            self.assertListEqual(self.profile[i], self.items[i][1])

    def test_columns(self):
        # Check if columns property returns correct values
        self.assertListEqual(self.profile.columns, self.col_names)

    def test_get(self):
        # Check if get method raises the correct exceptions
        message = 'identifier cannot be found in he Profile object.'
        with self.assertRaises(KeyError, msg=message):
            self.profile.get('f')
        with self.assertRaises(KeyError, msg=message):
            self.profile.get(1)
        with self.assertRaises(KeyError, msg=message):
            self.profile.get('A')
        # Check if get returns the correct values
        for idx, identifier in enumerate(self.row_names):
            self.assertListEqual(self.profile.get(identifier),
                                 self.items[idx][1])

    def test_samples(self):
        # Check if samples() returns a ExpressionProfile object
        self.assertIsInstance(self.profile.samples([0, 10, 2]),
                              ExpressionProfile)
        # Check if samples return correct values
        for idx, col in enumerate(self.col_names):
            self.assertListEqual(list(self.profile.samples([idx])),
                                 list(self.profile.samples([col],
                                                           is_index=False)))
            self.assertListEqual(self.profile.samples([idx]).columns,
                                 [col])
            self.assertListEqual(list(self.profile.samples([idx]).keys()),
                                 self.row_names)
        self.assertListEqual(list(self.profile.samples([0, 1]).columns),
                             ['A', 'B'])
        self.assertListEqual(list(self.profile.samples([0, 1]).keys()),
                             self.row_names)
        self.assertEqual(list(self.profile.samples(self.col_names,
                                                   is_index=False).items()),
                         self.items)
        col_indices = list(range(self.NUM_COLUMNS))
        self.assertEqual(list(self.profile.samples(col_indices).items()),
                         self.items)

    def test__str__(self):
        self.assertEqual(str(self.profile), str(self.data))

    def test_concat(self):
        profile1 = ExpressionProfile(self.data.loc[:, list('ABCDEF')])
        profile2 = ExpressionProfile(self.data.loc[:, list('GHIJKL')])
        profile3 = profile1.concat(profile2)
        self.assertListEqual(profile3.columns, list('ABCDEFGHIJKL'))
        self.assertListEqual(list(profile3.keys()), list(profile1.keys()))
        self.assertListEqual(list(profile3.keys()), list(profile2.keys()))
        pd.testing.assert_frame_equal(self.data, profile3.profile)


if __name__ == '__main__':
    unittest.main()
