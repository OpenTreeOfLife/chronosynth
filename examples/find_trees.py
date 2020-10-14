import chronosynth
from chronosynth import chronogram


studies = chronogram.find_trees()
prod = set(studies)

chronogram.print_endpoint()

chronogram.set_dev()

chronogram.print_endpoint()


studies = chronogram.find_trees()
dev = set(studies)

print("{} chrongrams in production\n".format(len(prod)))
print("{} chrongrams in dev\n".format(len(dev)))

print(prod.difference(dev))

print(dev.difference(prod))
