"""This mocule contains Exception objects.

"""
__all__ = ['InvalidGenesetFoldChangeFile', 'InvalidGenesetRepositoryFile',
            'InvalidExpressionFile', 'MismatchedProfileException',
            'InvalidNumberOfInputSamples', 'InvalidNumberOfOutputSamples',
            'InvalidValueforStandardDistributionOfNoise',
            'InvalidNumberOfGenes', 'CaseControlProfileMissMatch',
            'InvalidContrastException', 'NonExistingExpressionException']


class InvalidGenesetFoldChangeFile(Exception):
    def __init__(self, message):
        super().__init__(message)


class InvalidGenesetRepositoryFile(Exception):
    def __init__(self, message):
        super().__init__(message)


class InvalidExpressionFile(Exception):
    def __init__(self, message):
        super().__init__(message)


class MismatchedProfileException(Exception):
    def __init__(self, message):
        super().__init__(message)


class InvalidNumberOfInputSamples(Exception):
    def __init__(self, message):
        super().__init__(message)


class InvalidNumberOfOutputSamples(Exception):
    def __init__(self, message):
        super().__init__(message)


class InvalidValueforStandardDistributionOfNoise(Exception):
    def __init__(self, message):
        super().__init__(message)


class InvalidNumberOfGenes(Exception):
    def __init__(self, message):
        super().__init__(message)


class CaseControlProfileMissMatch(Exception):
    def __init__(self, message):
        super().__init__(message)


class InvalidContrastException(Exception):
    def __init__(self, message):
        super().__init__(message)


class NonExistingExpressionException(Exception):
    def __init__(self, message):
        super().__init__(message)
