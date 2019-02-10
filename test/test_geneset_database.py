from silver.geneset_database import GenesetDatabase
from collections import OrderedDict
import unittest


class TestGeneSetDatabase(unittest.TestCase):
    def setUp(self):
        genesets_items = [('Geneset1', ['gene4', 'gene5', 'gene2']),
                          ('Geneset2', ['gene1', 'gene3', 'gene2', 'gene4',
                                        'gene5'])]

        info_itmes = [('Geneset1', 'first geneset'),
                      ('Geneset2', 'second geneset')]

        self.databse = OrderedDict(genesets_items)
        self.info = OrderedDict(info_itmes)

        self.gsdb = GenesetDatabase(self.databse, self.info)

    def test__init__(self):
        db = GenesetDatabase(self.databse)
        self.assertDictEqual(db.database, self.databse)
        self.assertIsNone(db.info)

        self.assertDictEqual(self.gsdb.database, self.databse)
        self.assertDictEqual(self.gsdb.info, self.info)

    def test__getitem__(self):
        self.assertListEqual(self.gsdb['Geneset1'], self.databse['Geneset1'])
        self.assertListEqual(self.gsdb['Geneset2'], self.databse['Geneset2'])

    def test__iter__(self):
        self.assertListEqual(list(iter(self.gsdb)), list(iter(self.databse)))

    def test_keys(self):
        self.assertListEqual(list(self.gsdb.keys()), list(self.databse.keys()))

    def test_values(self):
        self.assertListEqual(list(self.gsdb.values()),
                             list(self.databse.values()))

    def test_items(self):
        self.assertListEqual(list(self.gsdb.items()),
                             list(self.databse.items()))

    def test_fromGMT(self):
        address = 'test/data/dummy_genesets.gmt'
        gsdb = GenesetDatabase.fromGMT(address)
        self.assertDictEqual(gsdb.database, self.databse)
        self.assertDictEqual(gsdb.info, self.info)

    def test_validate(self):
        background = ['gene1', 'gene2', 'gene3', 'gene4']
        self.gsdb.clean(background)
        self.assertListEqual(self.gsdb['Geneset1'], ['gene4', 'gene2'])
        self.assertListEqual(self.gsdb['Geneset2'], ['gene1', 'gene3', 'gene2',
                                                     'gene4'])
        self.gsdb.clean(background, min_size=3)
        message = 'gene set with a size smaller than min_size must be excluded.'
        with self.assertRaises(KeyError, msg=message):
            self.gsdb['Geneset1']
        self.assertListEqual(self.gsdb['Geneset2'], ['gene1', 'gene3', 'gene2',
                                                     'gene4'])
        self.gsdb.clean(background, max_size=3)
        message = 'gene set with a size greater than max_size must be excluded.'
        with self.assertRaises(KeyError, msg=message):
            self.gsdb['Geneset2']
