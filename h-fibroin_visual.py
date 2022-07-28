import re
import sys
from sys import argv
from math import ceil
from matplotlib import pyplot as plt, font_manager as fm
import numpy as np
import itertools

# For Plodia, uncomment lines 30-32, 66-67 (comment 64-65), 88-93

fprop = fm.FontProperties(family = 'sans serif')

def readFile(path):
	# Takes in path to fasta file with one or two h fibroin sequences

	with open(path) as file:
		sequences = file.read().strip().split("\n")
		seqs = []
		for seq in sequences:
			if not seq.startswith(">"):
				seqs.append(seq)
	return seqs

def organizeData(sequence):
	# Split given sequence at S blocks

	# print("yes")
	sequence = sequence.replace("-","")

	# # Caddisfly
	search = r"S.S.S"
	search = r"SISR" # hesp
	# search = r"QTPTI" # arcto
	# search = r"GPWGR"
	# search = r"KGRRG"

	# Spiders
	# MaSp2
	# search = r".GYGP"
	# AgSp2
	# search = r"[PS][EG][STA]TP"

	# Plodia
	# search = r"SA..A"

	# Trichoptera regex
	s_patterns = re.findall(search, sequence)
	i_patterns = re.split(search, sequence)

	
	first = i_patterns.pop(0)
	a_patterns = [first]
	for i in range(len(s_patterns)):
		pattern = s_patterns[i]  + i_patterns[i]
		a_patterns.append(pattern)

	return a_patterns

def gapPattterns(a_patterns):
	gap_patterns = []
	for pattern in a_patterns:
		if "-" in pattern:
			search = r"-+"
			split = re.split(search, pattern)
			findall = re.findall(search, pattern)
			# if len(findall[0]) < 10:
			# 	gap_patterns.append(pattern.replace("-",""))
			# 	continue
			if split[0] != "":
				gap_patterns.append(split[0])
			gap_patterns.append(findall[0])
			if split[1] != "":
				gap_patterns.append(split[1])
		else:
			gap_patterns.append(pattern)

	return gap_patterns

def removeGaps(gap_patterns):
	no_gaps = []
	for pattern in gap_patterns:
		if "-" not in pattern:
			no_gaps.append(pattern)
	return no_gaps

def gapIndex(gap_patterns, average):
	index = []
	subtract = 0
	for i, pattern in enumerate(gap_patterns):
		if "-" in pattern:
			rep = ceil(len(pattern) / average)
			for j in range(rep):
				index.append(i-subtract)
			subtract += 1
	return index

def addGaps(my_list, item, index):
	for i in index[::-1]:
		my_list.insert(i, item)
	return my_list

def organizeDataNumeric(a_patterns, number_dict, start):

	for pattern in a_patterns:
		if pattern not in number_dict:
			number_dict[pattern] = start 
			start += 1

	number_list = []

	for pattern in a_patterns:
		number_list.append(number_dict[pattern])

	return number_list, number_dict, start

def colorListSingle(a_patterns):
	# Build color dictionary and list for plotting

	possible_colors = ["#482677FF", "#FDE725FF", "#453781FF", "#fdc827",
	"#440154FF", "#55C667FF", "#DCE319FF", "#fad824", "#B8DE29FF",  "#0d0887", 
	"#238A8DFF", "#350498", "#33638DFF", "#95D840FF", "#2D708EFF", 
	"#73D055FF", "#404788FF", "#90d743", "#dae319", "#1F968BFF", "#3CBB75FF"]

	# # Alternate color palette - plodia
	# possible_colors = ["#fde725","#8b0aa5","#fada24","#e97158",
	# "#0d0887","#eb5760","#febd2a","#470b6a","#9f2a63",
	# "#b83289","#5302a3","#f0f921", "#d44842","#482878",
	# "#db5c68","#a31e9a","#350498","#f48849",
	# "#6f00a8","#fba238","#cc4778","#221150", "#f3e55d"]

	colors_key = ["#17202A"]
	hatches_key = [""]

	color_dict = {a_patterns[0]:"#17202A"}
	colors = ["#17202A"]
	
	hatches_dict = {a_patterns[0]:""}
	hatches = [""]

	i = 0
	hatch = ""

	x_patterns = [a_patterns[0]]

	for item in a_patterns[1:-1]:
		if item not in color_dict:
			similars = checkSimilarity(list(color_dict.keys()), item)
			if similars != None:
				similar = similars[0]
				new = similars[1]

				val = color_dict[similar]
				color_dict.pop(similar,None)
				color_dict[new] = val

				val2 = hatches_dict[similar]
				hatches_dict.pop(similar,None)
				hatches_dict[new] = val2

				colors.append(color_dict[new])
				hatches.append(hatches_dict[new])

				x_patterns.append(new)
				x_patterns = changeSimilars(new, x_patterns)

			else:
				x_patterns.append(item)
				color_dict[item] = possible_colors[i]
				colors.append(possible_colors[i])

				hatches_dict[item] = hatch
				hatches.append(hatch)

				colors_key.append(possible_colors[i])
				if colors_key.count(possible_colors[i]) > 1:
					hatches_key.append("..")
				else:
					hatches_key.append("")

				i += 1
				if i == len(possible_colors):
					i = 0
					hatch = ".."
		else:
			colors.append(color_dict[item])
			hatches.append(hatches_dict[item])
			x_patterns.append(item)

		
	x_patterns.append(a_patterns[-1])
	colors_key.append("#17202A")
	hatches_key.append("")

	color_dict[a_patterns[-1]] = "#17202A"
	colors.append("#17202A")
	hatches.append("")

	return x_patterns, colors_key, hatches_key, color_dict, colors, hatches

def colorListDouble(a_patterns1, a_patterns2):
	# Build two color lists using one combined dictionary

	possible_colors = ["#482677FF", "#FDE725FF", "#453781FF", "#fdc827",
	"#440154FF", "#55C667FF", "#DCE319FF", "#fad824", "#B8DE29FF",  "#0d0887", 
	"#238A8DFF", "#350498", "#33638DFF", "#95D840FF", "#2D708EFF", 
	"#73D055FF", "#404788FF", "#90d743", "#dae319", "#1F968BFF", "#3CBB75FF"]

	colors_key = ["#17202A"]
	hatches_key = [""]

	color_dict = {a_patterns1[0]:"#17202A", a_patterns2[0]:"#17202A"}
	colors1 = ["#17202A"]
	colors2 = ["#17202A"]
	
	hatches_dict = {a_patterns1[0]:"", a_patterns2[0]:""}
	hatches1 = [""]
	hatches2 = [""]

	i = 0
	hatch = ""

	x_patterns = [a_patterns1[0]]

	for item in a_patterns1[1:-1]:
		if item not in color_dict:
			similars = checkSimilarity(list(color_dict.keys()), item)
			if similars != None:
				similar = similars[0]
				new = similars[1]

				val = color_dict[similar]
				color_dict.pop(similar,None)
				color_dict[new] = val

				val2 = hatches_dict[similar]
				hatches_dict.pop(similar,None)
				hatches_dict[new] = val2

				colors1.append(color_dict[new])
				hatches1.append(hatches_dict[new])

				x_patterns.append(new)
				x_patterns = changeSimilars(new, x_patterns)

			else:
				x_patterns.append(item)
				color_dict[item] = possible_colors[i]
				colors1.append(possible_colors[i])

				hatches_dict[item] = hatch
				hatches1.append(hatch)

				colors_key.append(possible_colors[i])
				if colors_key.count(possible_colors[i]) > 1:
					hatches_key.append("..")
				else:
					hatches_key.append("")

				i += 1
				if i == len(possible_colors):
					i = 0
					hatch = ".."
		else:
			colors1.append(color_dict[item])
			hatches1.append(hatches_dict[item])
			x_patterns.append(item)

	for item in a_patterns2[1:-1]:
		if item not in color_dict:
			similars = checkSimilarity(list(color_dict.keys()), item)
			if similars != None:
				similar = similars[0]
				new = similars[1]

				val = color_dict[similar]
				color_dict.pop(similar,None)
				color_dict[new] = val

				val2 = hatches_dict[similar]
				hatches_dict.pop(similar,None)
				hatches_dict[new] = val2

				colors2.append(color_dict[new])
				hatches2.append(hatches_dict[new])

				x_patterns.append(new)
				x_patterns = changeSimilars(new, x_patterns)

			else:
				x_patterns.append(item)
				color_dict[item] = possible_colors[i]
				colors2.append(possible_colors[i])

				hatches_dict[item] = hatch
				hatches2.append(hatch)

				colors_key.append(possible_colors[i])
				if colors_key.count(possible_colors[i]) > 1:
					hatches_key.append("..")
				else:
					hatches_key.append("")

				i += 1
				if i == len(possible_colors):
					i = 0
					hatch = ".."
		else:
			colors2.append(color_dict[item])
			hatches2.append(hatches_dict[item])
			x_patterns.append(item)
		
	x_patterns.append(a_patterns2[-1])
	colors_key.append("#17202A")
	hatches_key.append("")

	color_dict[a_patterns1[-1]] = "#17202A"
	color_dict[a_patterns2[-1]] = "#17202A"
	colors1.append("#17202A")
	colors2.append("#17202A")
	hatches1.append("")
	hatches2.append("")

	return x_patterns, colors_key, hatches_key, color_dict, colors1, colors2, hatches1, hatches2

def checkSimilarity(keys, seq):

	return None # Use this if we want identical matches, no Xs
	# Check the number of differences for each key
	# If similar key is found, return key, else return none

	for key in keys:
		diffs = 0
		if len(key) == len(seq):
			for i in range(len(key)):
				if key[i] != seq[i]:
					diffs += 1
		# if diffs == 1:
		if diffs <= 5 and diffs >= 1:
			if combineSeqs(key,seq).count("X") <= 5:
				return [key, combineSeqs(key,seq)]
	return None

def combineSeqs(seq1, seq2):
	newSeq = []
	for i in range(len(seq1)):
		if seq1[i] != seq2[i]:
			newSeq.append("X")
		else:
			newSeq.append(seq1[i])
	return "".join(newSeq)

def unique(my_list):
	new_list = []
	for x in my_list:
		if x not in new_list:
			new_list.append(x)
	return new_list

def lengthList(a_patterns):
	# Generate a list of lengths of the patterns to use as height of bars in plot

	length = []
	for item in a_patterns:
		length.append(len(item))

	return length

def split(in_list):
	#Split list in half

	list1 = in_list[:int((len(a_patterns)+1)/2)]
	list2 = in_list[int((len(a_patterns)+1)/2):]
	return list1, list2

def horizontalBarPlot(a_patterns, colors, lengths, hatches, path):
	# Plot one allele, horizontal
	path = "horz_" + path
	first_half, second_half = split(a_patterns)
	first_colors, second_colors = split(colors)
	first_len, second_len = split(lengths)
	first_hatch, second_hatch = split(hatches)
	
	plt.rcParams['pdf.fonttype'] = 42
	fig, ax = plt.subplots(1, 2)
	
	plot1 = ax[0]
	plot2 = ax[1]
	plot1.set_title('', fontproperties=fprop)
	plot2.set_title('', fontproperties=fprop)
	
	fig.set_figheight(8)
	y_pos1 = np.arange(len(first_half))[::-1]
	y_pos2 = np.arange(len(second_half))[::-1]

	plot1.barh(y_pos1, first_len, align='center',
		color = first_colors, hatch = first_hatch)
	plot1.set_xlim([0,max(lengths)+2])
	plot1.set_yticks([])
	plot1.set_xlabel("Number of Residues")


	plot2.barh(y_pos2, second_len, align='center',
		color = second_colors, hatch = second_hatch)
	plot2.set_xlim([0,max(lengths)+2])
	plot2.set_yticks([])
	plot2.set_xlabel("Number of Residues")


	
	# plt.show()
	plt.savefig(path)

def verticalBarPlotSingle(a_patterns, colors, lengths, path, hatches, numbers):
	# Plot one allele

	plt.rcParams['pdf.fonttype'] = 42
	fig, ax = plt.subplots()
	
	fig.set_figheight(6)
	fig.set_figwidth(35)
	x_pos = np.arange(len(a_patterns))

	ax.bar(x_pos, lengths, align='center', color = colors, hatch = hatches)
	ax.set_xticks([])

	# When position is a positive number, the numbers are at the bottom of the bars (something around position = 2)
	# When position is a negative number, they are under the bars (something around position = -10)
	position = 2
	
	# Size options: 'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'
	size = "large"

	for i in range(len(x_pos)):
		ax.text(i,position,numbers[i], color = "red", ha = "center", fontsize = size)

	# ax.set_xlabel("Length of Allele", )
	# ax.set_ylabel("Number of Residues")
	ax.set_title('', fontproperties=fprop)

	# plt.show()
	plt.savefig(path)

def verticalBarPlotDouble(colors1, colors2, lengths1, lengths2, hatches1, hatches2, path):
	# Plot two alleles, colors coordinating

	plt.rcParams['pdf.fonttype'] = 42
	fig, ax = plt.subplots(2,1)

	plot1 = ax[0]
	plot2 = ax[1]

	plot1.set_title('', fontproperties=fprop)
	plot2.set_title('', fontproperties=fprop)
	
	fig.set_figheight(7)
	fig.set_figwidth(35)

	x_pos1 = np.arange(1, len(colors1)+1)
	x_pos2 = np.arange(1, len(colors2)+1)

	maxY = max([max(lengths1), max(lengths2)]) + 5
	# maxY = 200
	maxX = max([max(x_pos1), max(x_pos2)]) + 5

	plot1.bar(x_pos1, lengths1, align='center', color = colors1, hatch = hatches1)
	plot1.set_ylim([-5,maxY])
	plot1.set_xlim([-5,maxX])
	plot1.set_xticks(x_pos1)
	plot1.set_ylabel("Number of Residues")
	plot1.tick_params(labelrotation=90)

	plot2.bar(x_pos2, lengths2, align='center', color = colors2, hatch = hatches2)
	plot2.set_ylim([-5,maxY])
	plot2.set_xlim([-5,maxX])
	plot2.set_xticks(x_pos2)
	plot2.set_xlabel("Length of Allele")
	plot2.set_ylabel("Number of Residues")
	plot2.tick_params(labelrotation=90)

	# plt.show()
	plt.savefig(path)

def addCounts(patterns, x_patterns, colors, hatches, filter):
	patterns_counts = []
	colors_new = []
	hatches_new = []
	if filter:
		for i, pattern in enumerate(patterns):
			if i == 0 or i == (len(patterns) - 1):
				pattern = f"{pattern} ({2})"
				patterns_counts.append(pattern)
				colors_new.append(colors[i])
				hatches_new.append(hatches[i])
			elif x_patterns.count(pattern) < 3:
				continue
			else:
				pattern = f"{pattern} ({x_patterns.count(pattern)})"
				patterns_counts.append(pattern)
				colors_new.append(colors[i])
				hatches_new.append(hatches[i])
	else:
		for i, pattern in enumerate(patterns):
			if i == 0 or i == (len(patterns) - 1):
				pattern = f"{pattern} ({2})"
			else:
				pattern = f"{pattern} ({x_patterns.count(pattern)})"
			patterns_counts.append(pattern)
			colors_new.append(colors[i])
			hatches_new.append(hatches[i])
	return patterns_counts, colors_new, hatches_new    

def makeLegend(patterns, colors, hatches, path, numbers):

	fig, ax = plt.subplots()
	
	fig.set_figwidth(.5)
	fig.set_figheight(13)
	plt.rcParams['pdf.fonttype'] = 42

	x_pos = list(itertools.repeat(3, len(colors)))

	ax.barh(patterns[::-1], x_pos, align='center', color=colors[::-1], hatch = hatches[::-1])

	right = ax.spines["right"]
	right.set_visible(False)
	left = ax.spines["left"]
	left.set_visible(False)
	top = ax.spines["top"]
	top.set_visible(False)
	bottom = ax.spines["bottom"]
	bottom.set_visible(False)

	ax.set_title('', fontproperties=fprop)
	ax.set_xlim([0,4.5])
	ax.set_xticks([])
	ax.tick_params(axis='y', which='major', pad=15)
	ax.yaxis.tick_right()
	
	# Size options: 'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'
	size = "large"

	numbers = numbers[::-1]
	for i in range(len(x_pos)):
		ax.text(-3,i, numbers[i], color = "red", ha = "center", va = "center", fontsize = size)

	plt.savefig(path, bbox_inches="tight")

def replace(my_list, old, new):
	for i,j in enumerate(my_list):
		if j == old:
			my_list[i] = new
	return my_list

def changeSimilars(new, x_patterns):
	for i, pattern in enumerate(x_patterns):
		similars = checkSimilarity([pattern], new)
		if similars != None:
			similar = similars[0]
			new = similars[1]
			replace(x_patterns, similar, new)
	return x_patterns

def splitNucleotides(a_patterns, nuc_sequence):

	nuc_patterns = []
	start = 0
	stop = 0
	lengths = lengthList(a_patterns)
	for i in lengths:
		stop = start + (i * 3)
		seq = nuc_sequence[start:stop]
		nuc_patterns.append(seq)
		start = stop
	return nuc_patterns

def addNumbers(nuc_seqs, number_list):
	combined = []
	for i, nuc in enumerate(nuc_seqs):
		combo = f"{nuc} ({number_list[i]})"
		combined.append(combo)
	return combined

def addNumbersKey(patterns):
	new = []
	for i, pattern in enumerate(patterns):
		new_pattern = f"{i+1}: {pattern}"
		new.append(new_pattern)
	return new

def outputPatterns(a_patterns, path, sep):
	with open(path, "w") as file:
		for pattern in a_patterns:
			file.write(str(pattern) + sep)

def countsDict(a_patterns1, a_patterns2):
	counts = {}
	for a_patterns in [a_patterns1, a_patterns2]:
		for pattern in a_patterns:
			if pattern not in counts:
				counts[pattern] = 1
			else:
				counts[pattern] += 1
	return counts

def getMax(num, counts):
	keys = []
	values = []
	while len(keys) < num:
		max_val = 0
		max_key = ""
		for key, value in counts.items():
			if value > max_val and key not in keys:
				max_key = key
				max_val = value
		keys.append(max_key)
		values.append(max_val)
	return keys

def indexDictionary(num, counts, a_patterns1, a_patterns2):
	max_motifs = getMax(num, counts)
	index_dicts = []
	for i, motif in enumerate(max_motifs):
		index_dict = {}
		motif_num = i + 1
		copy_num = 1
		allele_num = 0
		for a_patterns in [a_patterns1, a_patterns2]:
			for j, pattern in enumerate(a_patterns):
				if pattern == motif:
					header = f">rep{allele_num+1}_{copy_num}"
					copy_num += 1
					index = [allele_num, j]
					index_dict[header] = index
			allele_num += 1
		index_dicts.append(index_dict)
	return index_dicts

def outputMotifs(nucs1, nucs2, index_dicts, path):
	nucs = [nucs1, nucs2]
	for i, fasta_dict in enumerate(index_dicts):
		file = path + f"_motif{i+1}.fasta"
		with open(file, "w") as fa:
			for header, index in sorted(fasta_dict.items()):
				fa.write(header + "\n")
				seq = nucs[index[0]][index[1]]
				fa.write(seq + "\n")

def outputGapsToFix(gaps, path):
	with open(path, "w") as file:
		for i in range(len(gaps[0])):
			line = ",".join([str(gaps[0][i]), gaps[1][i], gaps[2][i]])
			file.write(line + "\n")

def readCSV(input_csv):
	lengths = []
	colors = []
	hatches = []
	with open(input_csv) as csv:
		for line in csv:
			items = line.strip().split(",")
			lengths.append(int(items[0]))
			colors.append(items[1])
			hatches.append(items[2])
	return lengths, colors, hatches

def toString(numbers):
	strings = []
	for number in numbers:
		strings.append(str(number))
	return strings


if __name__ == "__main__":
	filename = argv[1]
	if filename.endswith(".csv"):
		if "lengend" in filename:
			keyCSV(input_csv)
		else:
			lengths1, colors1, hatches1 = readCSV(filename)
			lengths2, colors2, hatches2 = readCSV(argv[2])
			figure_path = filename.split(".")[0].strip("1") + ".pdf"
			verticalBarPlotDouble(colors1, colors2, lengths1, lengths2, hatches1, hatches2, figure_path)
	else:
		seqs = readFile(filename)
		if len(seqs) == 1:
			outfile = filename.split(".")[0] + "_patterns.txt"
			figure_path = filename.split(".")[0] + ".pdf"
			legend_path = filename.split(".")[0] + "_legend.pdf"
			seq = seqs[0]
			a_patterns = organizeData(seq)
			x_patterns, colors_key, hatches_key, color_dict, colors, hatches = colorListSingle(a_patterns)
			patterns_unique = unique(x_patterns)
			lengths = lengthList(a_patterns)

			numbers, number_dict, start = organizeDataNumeric(a_patterns, {}, 1)

			outputPatterns(addNumbers(a_patterns, numbers), outfile, "\n")
			# outputPatterns(x_patterns, "x_" + outfile, "\n")
			# horizontalBarPlot(a_patterns, colors, lengths, hatches, figure_path)
			print(len(a_patterns), len(numbers))
			verticalBarPlotSingle(a_patterns, colors, lengths, figure_path, hatches, toString(numbers))
			patterns_unique, colors_key, hatches_key = addCounts(patterns_unique, x_patterns, colors_key, hatches_key, False)
			makeLegend(patterns_unique, colors_key, hatches_key, legend_path, unique(numbers))

			legend_txt = filename.split(".")[0] + "_legend.txt"

			outputPatterns(addNumbersKey(patterns_unique), legend_txt, "\n")
			
		else:
			outfile1 = filename.split(".")[0] + "_patterns1.txt"
			outfile2 = filename.split(".")[0] + "_patterns2.txt"
			figure_path = filename.split(".")[0] + ".pdf"
			legend_path = filename.split(".")[0] + "_legend.pdf"
			seq1 = seqs[0]
			seq2 = seqs[1]
			a_patterns1 = organizeData(seq1.replace("-",""))
			a_patterns2 = organizeData(seq2.replace("-",""))

			x_patterns, colors_key, hatches_key, color_dict, colors1, colors2, hatches1, hatches2 = colorListDouble(a_patterns1, a_patterns2)
			lengths1 = lengthList(a_patterns1)
			lengths2 = lengthList(a_patterns2)
			patterns_unique = unique(x_patterns)

			numbers1, number_dict, start = organizeDataNumeric(a_patterns1, {}, 1)
			numbers2, number_dict, start = organizeDataNumeric(a_patterns2, number_dict, start)

			outputPatterns(addNumbers(a_patterns1, numbers1), outfile1, "\n")
			outputPatterns(addNumbers(a_patterns2, numbers2), outfile2, "\n")
			
			# outputPatterns(x_patterns, "x_" + outfile1, "\n")

			verticalBarPlotDouble(colors1, colors2, lengths1, lengths2, hatches1, hatches2, figure_path)
			patterns_unique, colors_key, hatches_key = addCounts(patterns_unique, x_patterns, colors_key, hatches_key, False)
			makeLegend(patterns_unique, colors_key, hatches_key, legend_path)

			legend_txt = filename.split(".")[0] + "_legend.txt"

			outputPatterns(addNumbersKey(patterns_unique), legend_txt, "\n")

			fixPath1 = filename.split(".")[0] + "_fix_gaps_1.csv"
			fixPath2 = filename.split(".")[0] + "_fix_gaps_2.csv"
			# outputGapsToFix([lengths1, colors1, hatches1],fixPath1)
			# outputGapsToFix([lengths2, colors2, hatches2],fixPath2)
			
			
