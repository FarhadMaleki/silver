import unittest
import numpy as np
import pandas as pd
from silver.repository import Repository
from silver.dexpress import TTestDExpress
from silver.expression_profile import ExpressionProfile
from silver.exceptions import NonExistingExpressionException


class TestTTestDExpress(unittest.TestCase):
    def setUp(self):
        np.random.seed(213456)
        # Read the expression dataset
        address = 'test/data/dummy_expression.txt'
        data = ExpressionProfile(pd.read_csv(address, sep='\t', index_col='ID'))
        # Build control and case expression profiles
        num_ctrls = 6
        num_cases = 6
        ctrl_indices = list(range(num_ctrls))
        case_indices = list(range(num_ctrls, num_ctrls + num_cases))
        self.ctrl_profile = data.samples(ctrl_indices)
        self.case_profile = data.samples(case_indices)

        # Create a repository
        num_simulated_cases = 3
        self.repository = Repository(self.case_profile,
                                     num_sim_cases=num_simulated_cases,
                                     num_repetitions=1,
                                     random_state=123456)
        alpha = 0.05

        self.ttest_dexpress = TTestDExpress(self.repository, alpha)

    def test_log_fold_change_of_1(self):
        ctrl_expression = [1.12, 1.13, 0.97]
        lower, upper = 0.5, 1.8
        case_expression = self.ttest_dexpress(ctrl_expression, (lower, upper))
        correct = set(case_expression).issubset(set(self.case_profile.get('b')))
        self.assertTrue(correct)

    def test_log_fold_change_of_minus_1(self):
        ctrl_expression = [1.95, 1.97, 2.18]
        lower, upper = -1.5, -0.5
        case_expression = self.ttest_dexpress(ctrl_expression, (lower, upper))
        correct = set(case_expression).issubset(set(self.case_profile.get('a')))
        self.assertTrue(correct)

    def test_log_fold_change_of_minus_2(self):
        ctrl_expression = [3.05, 3.17, 3.38]
        lower, upper = -2.5, -1.5
        case_expression = self.ttest_dexpress(ctrl_expression, (lower, upper))
        correct = set(case_expression).issubset(set(self.case_profile.get('a')))
        self.assertTrue(correct)

    def test_raising_NonExistingExpressionException(self):
        msg = 'Cannot find expression measures with requested fold changes.'
        # Check if NonExistingExpressionException is raised when it is expected
        with self.assertRaises(NonExistingExpressionException, msg=msg):
            ctrl_expression = [1.05, 1.17, 0.38]
            lower, upper = -2.5, -1.5
            self.ttest_dexpress(ctrl_expression, (lower, upper), force=False)
        # Check if NonExistingExpressionException is raised when it is expected
        with self.assertRaises(NonExistingExpressionException, msg=msg):
            ctrl_expression = [5., 4.87, 4.98]
            lower, upper = 1.0, 1.5
            self.ttest_dexpress(ctrl_expression, (lower, upper), force=False)

    def test_force_truein_ttest_dexpress_call(self):
        ctrl_expression = [1.05, 1.17, 0.38]
        lower, upper = -2.5, -1.5
        expression = self.ttest_dexpress(ctrl_expression,
                                         (lower, upper),
                                         force=True)
        self.assertTrue(np.mean(expression) >= np.mean(ctrl_expression) + lower)
        self.assertTrue(np.mean(expression) <= np.mean(ctrl_expression) + upper)



if __name__ == '__main__':
    unittest.main()
