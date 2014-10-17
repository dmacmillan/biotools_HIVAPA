#!/home4/dmacmill/python/Python-2.7.2/python
import cgi, sys, math
from decimal import Decimal

sys.stderr = open("/home4/dmacmill/public_html/oakdrum/HIVAPA/error-cgi.log", "a")

form = cgi.FieldStorage()
hlas = form.getvalue("hlasname")
seqs = form.getvalue("sequencesname")
runHIVAPA = form.getvalue("runHIVAPA")

def printHtmlHeaders():
	print "Content-Type: text/html"
	print
	print """<!DOCTYPE html><html><head>
	<script src="//ajax.googleapis.com/ajax/libs/jquery/1.11.1/jquery.min.js"></script>
	<script src="../cgi-bin/script.js"></script>
	<link rel="stylesheet" href="../css/style.css"></head><body>"""

def printFileHeaders(filename):
	print "Content-Disposition: attachment; filename=\""+filename+"\""
	print "Content-Type:application/octet-stream; name=\""+filename+"\""
	print

codon_dict = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
			  'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
			  'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
			  'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W',
			  'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
			  'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
			  'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
			  'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
			  'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
			  'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
			  'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
			  'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
			  'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
			  'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
			  'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
			  'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
			  '---':'-', 'XXX':'-', '???':'?'}

mixture_dict = {'W':'AT', 'R':'AG', 'K':'GT', 'Y':'CT',
				'S':'CG', 'M':'AC', 'V':'AGC', 'H':'ATC',
				'D':'ATG', 'B':'TGC', 'N':'ATGC', '-':'-'}

def resolveCodon(codon):
	nonmix = []
	if (codon in codon_dict):
		return [codon]
	elif (codon.count('-') + codon.count('X') == 3):
		return ['---']
	elif (1 <= codon.count('-') <= 2) or (1 <= codon.count('X') <= 2):
		return ['???']
	for base in codon:
		# Check for mixtures
		if (base in mixture_dict):
			if (not nonmix):
				nonmix = [x for x in mixture_dict[base]]
			else:
				nonmix = [x+y for x in nonmix for y in mixture_dict[base]]
		else:
			if (not nonmix):
				nonmix.append(base)
			else:
				nonmix = [x+base for x in nonmix]
	return nonmix

# Flag can be 0, 1, or 2 depending on the desired output
# Flag = 1 will output all mixtures as "X"
# Flag = 2 will output all synonymous mixtures as they are and all non-synonymous mixtures as "X"
# Flag = 3 will output all mixtures in the format [A/B] if a mixture encodes for amino acid A or B
def translateDNA(sequence, resolvecharacter="X", flag=2):
	sequence = sequence.translate(None, ' \n\r\n').upper()
	aaseq = []
	i = 0
	while i < len(sequence):
		codon = resolveCodon(sequence[i:i+3])
		# If the codon has no mixture bases just add it to the amino acid chain
		if len(codon) <= 1:
			aaseq.append(codon_dict[codon[0]])
		# Codon contains mixture base
		else:
			# If flag is set to 1
			if (flag == 1):
				aaseq.append(resolvecharacter)
			# If flag is set to 2
			elif (flag == 2):
				unique = set([codon_dict[potential] for potential in codon])
				# If there is more than resolved one amino acid
				if (len(unique) > 1):
					aaseq.append(resolvecharacter)
				else:
					aaseq.append(unique.pop())
			# If flag is set to 3
			else:
				unique = set([codon_dict[potential] for potential in codon])
				# If there is more than resolved one amino acid
				if (len(unique) > 1):
					aaseq.append('['+('/').join(unique)+']')
				else:
					aaseq.append(unique.pop())
		i += 3
	return aaseq

def parse(inputText):
	val = [x.split('\t') for x in inputText.splitlines()]
	return val

def parseHLA(hla, res=4):
	rval = hla.strip()
	rval = hla.translate(None, ":*")
	try:
		int(rval[-1])
	except (ValueError, IndexError) as e:
		rval = rval[:-1]
	return rval[:res+1]

def groupHLA(hlas):
	rdic = {}
	for pair in hlas:
		hla = parseHLA(pair[0])
		loc = pair[1][:-1]
		aa = pair[1][-1]
		if hla not in rdic:
			rdic[hla] = [[loc, aa]]
		else:
			rdic[hla].append([loc, aa])
	return rdic

def getSeqs(seqs):
	d = {}
	for patient in seqs:
		pid = patient[0]
		nuseq = patient[-1]
		if (len(nuseq) % 3 != 0):
			aaseq = "Not divisible by 3"
		else:
			aaseq = translateDNA(nuseq)
		d[pid] = {'A': [], 'B': [], 'C': [], 'seq': aaseq}
		for hla in patient[1:-1]:
			hla = parseHLA(hla)
			if (hla == ""):
				continue
			if (hla[0].upper() == 'A'):
				d[pid]['A'].append(hla)
			elif (hla[0].upper() == 'B'):
				d[pid]['B'].append(hla)
			else:
				d[pid]['C'].append(hla)
	return d

def getResults(patients, groupedHLAs):
    results = {}
    for patient in patients:
        potential, actual = 0, 0
        results[patient] = [{},{}]
        #N = len(patients)
        for hla in groupedHLAs:
            patientHLAs = patients[patient][hla[0].upper()]
            if (hla[:len(hla)] in [x[:len(hla)] for x in patientHLAs]):
                potential = len(groupedHLAs[hla])
                results[patient][0][hla] = groupedHLAs[hla]
            elif any(x[:len(x)] == hla[:len(x)] for x in patientHLAs) or (len(patientHLAs) < 2):
                #N -= 1
                continue
            for mut in groupedHLAs[hla]:
                try:
                    patientAA = patients[patient]['seq'][int(mut[0])-1]
                except IndexError:
                    #N -= 1
                    continue
                if (len(patientAA) > 1):
                    patientAA = patientAA[1:-1].split('/')
                else:
                    patientAA = [patientAA]
                if (hla in patientHLAs) and (mut[1] in patientAA):
                    actual += 1
                    if (hla not in results[patient][1]):
                        results[patient][1][hla] = [mut]
                    else:
                        results[patient][1][hla].append(mut)
        #results[patient].append(N)
    return results  

def computeMedian(results):
    N = len(results)
    r = {}
    for patient in results:
        r[patient] = []
        for patient2 in results:
            actual, potential = 0, 0
            if (patient == patient2):
                N -= 1
                continue
            for hla in results[patient][0]:
                if (hla in results[patient2][0]):
                    potential += len(results[patient][0][hla])
                for mut in results[patient][0][hla]:
                    if (hla in results[patient2][1]) and (mut in results[patient2][1][hla]):
                        actual += 1
            r[patient].append([actual, potential])
        r[patient] = sorted(r[patient], key = lambda(x): (float(x[0])/x[1]))
    return r
            

def displayResults(results):
    medians = computeMedian(results)
    print '''<table id="output_table">
             <th>patientID</th>
             <th>potential</th>
             <th>actual</th>
             <th>N</th>
             <th>median</th>
             <th>IQR</th>
             <th>max</th>
             <th>min</th>'''
    for patient in results:
        print '<tr>'
        print '<td>{}</td>'.format(patient)
        potential = 0
        for hla in results[patient][0]:
            potential += len(results[patient][0][hla])
        print '<td>{}</td>'.format(potential)
        print '<td>'
        for hla in results[patient][1]:
            print '{}: ({})<br>'.format(hla, (',').join([('').join(x) for x in results[patient][1][hla]]))
        print '</td>'
        print '<td>N (to be implemented)</td>'
        if (len(medians[patient]) % 2 != 0):
            median = medians[patient][len(medians[patient])/2]
        else:
            upper = medians[patient][len(medians[patient])/2]
            upper = float(upper[0]) / upper[1]
            lower = medians[patient][(len(medians[patient])/2)-1]
            lower = float(lower[0]) / lower[1]
            median = (upper + lower) / 2
        print '<td>{}</td>'.format(median)
        q1 = medians[patient][len(medians[patient])/4]
        q3 = medians[patient][(len(medians[patient])/4)*3]
        IQR = (float(q3[0])/q3[1]) - (float(q1[0])/q1[1])
        print '<td>{}</td>'.format(IQR)
        print '<td>{}</td>'.format(medians[patient][-1])
        print '<td>{}</td>'.format(medians[patient][0])
    print '</table>'

if (runHIVAPA is not None):
    printHtmlHeaders()
    hlas = parse(hlas)
    groupedHLAs = groupHLA(hlas)
    seqs = parse(seqs)
    patients = getSeqs(seqs)
    results = getResults(patients, groupedHLAs)
    displayResults(results)
