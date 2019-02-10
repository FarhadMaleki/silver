from collections import OrderedDict


class GenesetDatabase(object):
    def __init__(self, database, info=None):
        assert isinstance(database, OrderedDict)
        self.database = database
        self.info = info

    def __getitem__(self, identifier):
        return self.database[identifier]

    def __iter__(self):
        return iter(self.database)

    def keys(self):
        return self.database.keys()

    def items(self):
        return self.database.items()

    def values(self):
        return self.database.values()

    @classmethod
    def fromGMT(cls, address):
        sep = '\t'
        database, info = OrderedDict(), OrderedDict()
        with open(address) as fin:
            for line in fin:
                line = line.strip()
                if line == '':
                    continue
                words = [w.strip() for w in line.split(sep)]
                name = words[0]
                info[name] = words[1]
                database[name] = words[2:]
        return GenesetDatabase(database, info)

    def clean(self, background, min_size=1, max_size=float('inf')):
        background = set(background)
        db = OrderedDict()
        for name, gs in self.database.items():
            geneset = [g for g in gs if g in background]
            if min_size <= len(geneset) <= max_size:
                db[name] = geneset
            else:
                del self.info[name]
        self.database = db

