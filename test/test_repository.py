import unittest
import pandas as pd
from silver.repository import Repository
from silver.expression_profile import ExpressionProfile


class TestRepository(unittest.TestCase):
    def setUp(self):
        address = 'test/data/dummy_expression.txt'
        data = ExpressionProfile(pd.read_csv(address, index_col='ID', sep='\t'))
        self.data_length = len(data)
        self.num_ctrls, self.num_cases = 6, 6
        case_indices = list(range(self.num_ctrls,
                                  self.num_ctrls+self.num_cases))
        self.case_profile = data.samples(case_indices)
        self.num_sim_cases = 3
        self.repository_1rep = Repository(self.case_profile,
                                          num_sim_cases=self.num_sim_cases,
                                          num_repetitions=1)
        self.repository_3rep = Repository(self.case_profile,
                                          num_sim_cases=self.num_sim_cases,
                                          num_repetitions=3)

    def test__len__(self):
        self.assertEqual(len(self.repository_1rep), self.data_length)
        self.assertEqual(len(self.repository_3rep), self.data_length* 3)

    def test__getitem__(self):
        for i in range(self.data_length):
            self.assertEqual(len(self.repository_1rep[i]), self.num_sim_cases)
            self.assertEqual(set(self.repository_1rep[i]) -
                             set(self.case_profile[i]), set())
        with self.assertRaises(IndexError, msg='Index out of bound.'):
            self.repository_1rep[len(self.repository_1rep)]

        with self.assertRaises(IndexError, msg='Index out of bound.'):
            self.repository_1rep[-1]

        with self.assertRaises(IndexError, msg='Index out of bound.'):
            self.repository_3rep[len(self.repository_3rep)]

        with self.assertRaises(IndexError, msg='Index out of bound.'):
            self.repository_3rep[-1]

    def test__iter__(self):
        i = 0
        for val in self.repository_1rep:
            self.assertEqual(len(val), self.num_sim_cases)
            self.assertEqual(set(val) - set(self.case_profile[i]), set())
            i += 1
        i = 0
        for idx, val in enumerate(self.repository_3rep):
            i = idx % self.data_length
            self.assertEqual(len(val), self.num_sim_cases)
            self.assertEqual(set(val) - set(self.case_profile[i]), set())


if __name__ == '__main__':
    unittest.main()
