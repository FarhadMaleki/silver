import unittest
import numpy as np
import pandas as pd
from silver.expression_profile import ExpressionProfile
from silver.dataset import Dataset
from silver.repository import Repository
from silver.dexpress import TTestDExpress


class TestDataset(unittest.TestCase):
    def setUp(self):
        # Read the expression dataset
        address = 'test/data/dummy_expression.txt'
        data = pd.read_csv(address, sep='\t', index_col='ID')
        # Build control and case expression profiles
        self.num_ctrls = 6
        self.num_cases = 6
        ctrl_indices = list(range(self.num_ctrls))
        case_indices = list(range(self.num_ctrls, (self.num_ctrls +
                                                   self.num_cases)))
        self.ctrl_profile = ExpressionProfile(data.iloc[:, ctrl_indices])
        self.case_profile = ExpressionProfile(data.iloc[:, case_indices])
        # Create a dataset

    def test_make_custom_replicates(self):
        num_sim_ctrls, num_sim_cases = 3, 3
        sim_ctrl_indices = list(range(num_sim_ctrls))
        sim_case_indices = list(range(num_sim_ctrls, self.num_ctrls))
        sim_ctrl_columns = [self.ctrl_profile.columns[i]
                            for i in sim_ctrl_indices]
        sim_case_columns = [self.ctrl_profile.columns[i]
                            for i in sim_case_indices]

        dataset = Dataset(self.ctrl_profile, self.case_profile)

        ctrls, cases = dataset.make_custom_replicates(sim_ctrl_indices,
                                                      sim_case_indices)
        self.assertListEqual(ctrls.columns, sim_ctrl_columns)
        self.assertListEqual(cases.columns, sim_case_columns)
        self.assertListEqual(list(ctrls.keys()), list(cases.keys()))
        self.assertEqual(len(ctrls.columns), num_sim_ctrls)
        for i, val in enumerate(ctrls):
            self.assertTrue(set(val).issubset(set(self.ctrl_profile[i])))
            self.assertEqual(len(val), num_sim_ctrls)

        self.assertEqual(len(cases.columns), num_sim_cases)
        for i, val in enumerate(cases):
            self.assertTrue(set(val).issubset(set(self.ctrl_profile[i])))
            self.assertEqual(len(val), num_sim_cases)

    def test_diff_express(self):
        # Instantiate the Dataset object
        num_sim_ctrls, num_sim_cases = 3, 3
        sim_ctrl_indices = list(range(num_sim_ctrls))
        sim_case_indices = list(range(num_sim_ctrls, self.num_ctrls))
        # sim_ctrl_columns = [self.ctrl_profile.columns[i]
        #                     for i in sim_ctrl_indices]
        # sim_case_columns = [self.ctrl_profile.columns[i]
        #                     for i in sim_case_indices]

        dataset = Dataset(self.ctrl_profile, self.case_profile)
        ctrls, cases = dataset.make_custom_replicates(sim_ctrl_indices,
                                                      sim_case_indices)
        # Create a repository
        num_sim_cases = 3
        repository = Repository(self.case_profile,
                                num_sim_cases=num_sim_cases,
                                num_repetitions=2)
        # instantiate TTestDExpress object and apply differential expression
        ttest_express = TTestDExpress(repository, alpha=0.05)
        fc = {'a': (.5, 1.5)}
        sim_ctrls, sim_cases = dataset.diff_express(ctrls, cases, fc,
                                                    ttest_express)
        self.assertIs(ctrls, sim_ctrls)
        self.assertListEqual(cases.columns, sim_cases.columns)
        self.assertListEqual(list(cases.keys()), list(sim_cases.keys()))
        for gene in set(set(cases.keys() - fc.keys())):
            difference = sum(np.abs(np.array(cases.get(gene)) -
                                    np.array(sim_cases.get(gene))))
            self.assertAlmostEqual(difference, 0, delta=10E-7)
        for gene in set(fc.keys()):
            ctrl_median = np.median(sim_ctrls.get(gene))
            case_median = np.median(sim_cases.get(gene))
            lower_bound = fc[gene][0] + ctrl_median
            upper_bound = fc[gene][1] + ctrl_median
            self.assertTrue(lower_bound <= case_median <= upper_bound)


if __name__ == '__main__':
    unittest.main()
