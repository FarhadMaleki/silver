import unittest
from silver.exceptions import *


class TestExceptions(unittest.TestCase):

    def test_InvalidGenesetFoldChangeFile(self):
        msg = 'An exception was raised.'
        try:
            raise InvalidGenesetFoldChangeFile(msg)
            self.assertTrue(False)
        except InvalidGenesetFoldChangeFile as e:
            self.assertTrue(str(e), msg)
        except:
            self.assertTrue(False)



    def test_InvalidGenesetRepositoryFile(self):
        msg = 'An exception was raised.'
        try:
            raise InvalidGenesetRepositoryFile(msg)
            self.assertTrue(False)
        except InvalidGenesetRepositoryFile as e:
            self.assertTrue(str(e), msg)
        except:
            self.assertTrue(False)



    def test_InvalidExpressionFile(self):
        msg = 'An exception was raised.'
        try:
            raise InvalidExpressionFile(msg)
            self.assertTrue(False)
        except InvalidExpressionFile as e:
            self.assertTrue(str(e), msg)
        except:
            self.assertTrue(False)



    def test_MismatchedProfileException(self):
        msg = 'An exception was raised.'
        try:
            raise MismatchedProfileException(msg)
            self.assertTrue(False)
        except MismatchedProfileException as e:
            self.assertTrue(str(e), msg)
        except:
            self.assertTrue(False)



    def test_InvalidNumberOfInputSamples(self):
        msg = 'An exception was raised.'
        try:
            raise InvalidNumberOfInputSamples(msg)
            self.assertTrue(False)
        except InvalidNumberOfInputSamples as e:
            self.assertTrue(str(e), msg)
        except:
            self.assertTrue(False)



    def test_InvalidNumberOfOutputSamples(self):
        msg = 'An exception was raised.'
        try:
            raise InvalidNumberOfOutputSamples(msg)
            self.assertTrue(False)
        except InvalidNumberOfOutputSamples as e:
            self.assertTrue(str(e), msg)
        except:
            self.assertTrue(False)



    def test_InvalidValueforStandardDistributionOfNoise(self):
        msg = 'An exception was raised.'
        try:
            raise InvalidValueforStandardDistributionOfNoise(msg)
            self.assertTrue(False)
        except InvalidValueforStandardDistributionOfNoise as e:
            self.assertTrue(str(e), msg)
        except:
            self.assertTrue(False)



    def test_InvalidNumberOfGenes(self):
        msg = 'An exception was raised.'
        try:
            raise InvalidNumberOfGenes(msg)
            self.assertTrue(False)
        except InvalidNumberOfGenes as e:
            self.assertTrue(str(e), msg)
        except:
            self.assertTrue(False)



    def test_CaseControlProfileMissMatch(self):
        msg = 'An exception was raised.'
        try:
            raise CaseControlProfileMissMatch(msg)
            self.assertTrue(False)
        except CaseControlProfileMissMatch as e:
            self.assertTrue(str(e), msg)
        except:
            self.assertTrue(False)



    def test_InvalidContrastException(self):
        msg = 'An exception was raised.'
        try:
            raise InvalidContrastException(msg)
            self.assertTrue(False)
        except InvalidContrastException as e:
            self.assertTrue(str(e), msg)
        except:
            self.assertTrue(False)



    def test_NonExistingExpressionException(self):
        msg = 'An exception was raised.'
        try:
            raise NonExistingExpressionException(msg)
            self.assertTrue(False)
        except NonExistingExpressionException as e:
            self.assertTrue(str(e), msg)
        except:
            self.assertTrue(False)
