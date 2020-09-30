# TEST manually getting conflict data for all studies with branchLengthMode=="time"


# chronogram.set_dev()
# chronogram.set_prod()


studies = chronogram.find_trees()

f_errors = open("examples/conflict_errors.txt", "w+")
f_no_conflict = open("examples/conflict_none.txt", "w+")
f_good = open("examples/coflict_good.txt", "w+")
for tag in studies:
    try:
        res = chronogram.map_conflict_ages(tag)
    except ValueError:
        f_errors.write(tag)
        f_errors.write('\n')
        print("study", tag, "threw an ERROR on map_conflict ages")
        continue
    if res==None:
        print("study", tag, "has no conflict data")
        f_no_conflict.write(tag)
        f_no_conflict.write('\n')
    else:
        print("study", tag, "has", len(res["supported_nodes"]), "supported_nodes")
        f_good.write(tag)
        f_good.write('\n')

f_errors.close()
f_no_conflict.close()
f_good.close()
