from chronosynth import chronogram

#sc = chronogram.find_trees(value='ot:substitutionCount')
#un = chronogram.find_trees(value='ot:undefined')
#cc = chronogram.find_trees(value='ot:changesCount')
#ot = chronogram.find_trees(value='ot:other')

source = 'ot_303@Tr62284'

maps = chronogram.map_conflict_nodes(source)


##Next create priors file using fastdate