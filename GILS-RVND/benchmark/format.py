filepath = "target_data"
with open(filepath) as fp:
	lines = fp.read().splitlines()
with open(filepath, "w") as fp:
	for line in lines:
		data = line.split(' : ')
		new_path = "instances/" + data[0] + ".tsp:" + data[1]
		print(new_path, file=fp)
