#!/usr/bin/env python3
import sys
fastaseq = open(sys.argv[1], "r")
genes = {}

for line in fastaseq:
	line = line.rstrip()
	if line[0]=='>':
		spacepos = line.find(' ')
		gene = line[:spacepos].strip('>')
		if gene not in genes.keys():
			genes[gene] = ''
	else:
		genes[gene] += line.strip('\n')

# para determinar a sequência reversa complementar de cada sequência, e fazer as traduções para aminoácidos, podemos usar dicionários
complementBase = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
translation_table = {
    'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
    'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R',
    'AAT':'N', 'AAC':'N',
    'GAT':'D', 'GAC':'D',
    'TGT':'C', 'TGC':'C',
    'CAA':'Q', 'CAG':'Q',
    'GAA':'E', 'GAG':'E',
    'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
    'CAT':'H', 'CAC':'H',
    'ATT':'I', 'ATC':'I', 'ATA':'I',
    'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
    'AAA':'K', 'AAG':'K',
    'ATG':'M',
    'TTT':'F', 'TTC':'F',
    'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
    'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S',
    'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    'TGG':'W',
    'TAT':'Y', 'TAC':'Y',
    'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
    'TAA':'*', 'TGA':'*', 'TAG':'*'
}

# podemos criar uma lista para ser usada posteriormente para determinar a ORF mais longa
longest_orf = []

for gene in genes.keys():
# como já feito em um exercício anterior, podemos criar listas para armazenar os códons em forward, a sequência reversa e sua complementar
# o bloco em seguida irá definir a sequência reversa, estabelecer a fita complementar usando o dicionário
# e finalmente, criar a lista da fita complementar com a sequência separada em códons
	codons1 = []
	reverseSeq = []
	complementSeq = []
	complementar = []
	for i in range(0, len(genes[gene]), 3):
		codons1.append(genes[gene][i:i+3])	
	reverseSeq.append(genes[gene][::-1])
	reverseSeq = ''.join(reverseSeq)
	for base in reverseSeq:
		complementSeq.append(complementBase[base])
	complementSeq = ''.join(complementSeq)
	for i in range (0, len(complementSeq), 3):
		complementar.append(complementSeq[i:i+3])
# em seguida, encontraremos o códon de iniciação 'ATG' nas sequências
# depois, encontramos algum dos códons de parada 'TAA', 'TAG' ou 'TGA' para definir a ORF
	longest_orf = ['0']
	orf1 = []
	if 'ATG' in codons1:
		start1 = codons1. index('ATG')
		if 'TAA' in codons1:
			end1 = codons1.index('TAA')
		elif 'TGA' in codons1:
			end1 = codons1.index('TGA')
		elif 'TAG' in codons1:
			end1 = codons1.index('TAG')
		else:
			end1 = len(codons1)-1
		for i in range (start1, end1+1):
			orf1.append(codons1[i])
# vamos considerar que a primeira ORF é a mais longa
		orf1 = ''.join(orf1)
		longest_orf = []
		longest_orf .append(orf1)
		longest_orf.append(gene+'_frame1_'+str(start1)+'_'+str(end1))
		
# repetimos o processo para a fita reversa
	orf1rev = []
	if 'ATG' in complementar:
		start1rev = complementar. index('ATG')
		if 'TAA' in complementar:
			end1rev = complementar.index('TAA')
		elif 'TGA' in complementar:
			end1rev = complementar.index('TGA')
		elif 'TAG' in complementar:
			end1rev = complementar.index('TAG')
		else:
			end1rev = len(complementar)-1
		for i in range (start1rev, end1rev+1):
			orf1rev.append(complementar[i])
		orf1rev = ''.join(orf1rev)

		#caso a nova ORF seja maior que a anterior, vamos substituir na lista
		if len(orf1rev) > len(longest_orf[0]):
			longest_orf = []
			longest_orf .append(orf1rev)
			longest_orf.append(gene+'_frame2_'+str(start1rev)+'_'+str(end1rev))
		
#repetindo todo o processo para o segundo reading frame
	codons2 = []
	complementar2 = []
	for i in range(1, len(genes[gene]), 3):
		codons2.append(genes[gene][i:i+3])
	for i in range (1, len(complementSeq), 3):
		complementar2.append(complementSeq[i:i+3])	
	orf2 = []
	if 'ATG' in codons2:
		start2 = codons2. index('ATG')
		if 'TAA' in codons2:
			end2 = codons2.index('TAA')
		elif 'TGA' in codons2:
			end2 = codons2.index('TGA')
		elif 'TAG' in codons2:
			end2 = codons2.index('TAG')
		else:
			end2 = len(codons2)-1
		for i in range (start2, end2+1):
			orf2.append(codons2[i])
		orf2 = ''.join(orf2)
		if len(orf2) > len(longest_orf[0]):
			longest_orf = []
			longest_orf .append(orf2)
			longest_orf.append(gene+'_frame3_'+str(start2)+'_'+str(end2))
	
	orf2rev = []
	if 'ATG' in complementar2:
		start2rev = complementar2. index('ATG')
		if 'TAA' in complementar2:
			end2rev = complementar2.index('TAA')
		elif 'TGA' in complementar2:
			end2rev = complementar2.index('TGA')
		elif 'TAG' in complementar2:
			end2rev = complementar2.index('TAG')
		else:
			end2rev = len(complementar2)-1
		for i in range (start2rev, end2rev+1):
			orf2rev.append(complementar2[i])
		orf2rev = ''.join(orf2rev)
		if len(orf2rev) > len(longest_orf[0]):
			longest_orf = []
			longest_orf .append(orf2rev)
			longest_orf.append(gene+'_frame4_'+str(start2rev)+'_'+str(end2rev))

#repetindo o processo para o terceiro reading frame
	codons3 = []
	complementar3 = []
	for i in range(2, len(genes[gene]), 3):
		codons3.append(genes[gene][i:i+3])
	for i in range (2, len(complementSeq), 3):
		complementar3.append(complementSeq[i:i+3])
	
	orf3 = []
	if 'ATG' in codons3:
		start3 = codons3. index('ATG')
		if 'TAA' in codons3:
			end3 = codons3.index('TAA')
		elif 'TGA' in codons3:
			end3 = codons3.index('TGA')
		elif 'TAG' in codons3:
			end3 = codons3.index('TAG')
		else:
			end3 = len(codons3)-1
		for i in range (start3, end3+1):
			orf3.append(codons3[i])
		orf3 = ''.join(orf3)
		if len(orf3) > len(longest_orf[0]):
			longest_orf = []
			longest_orf .append(orf3)
			longest_orf.append(gene+'_frame5_'+str(start3)+'_'+str(end3))

	orf3rev = []
	if 'ATG' in complementar3:
		start3rev = complementar3. index('ATG')
		if 'TAA' in complementar3:
			end3rev = complementar3.index('TAA')
		elif 'TGA' in complementar3:
			end3rev = complementar3.index('TGA')
		elif 'TAG' in complementar3:
			end3rev = complementar3.index('TAG')
		else:
			end3rev = len(complementar3)-1
		for i in range (start3rev, end3rev+1):
			orf3rev.append(complementar3[i])
		orf3rev = ''.join(orf3rev)
		if len(orf3rev) > len(longest_orf[0]):
			longest_orf = []
			longest_orf .append(orf3rev)
			longest_orf.append(gene+'_frame6_'+str(start3rev)+'_'+str(end3rev))
			
	orf = longest_orf[0]
	bases = []
	pep = []
# separando a ORF em códons, passando para uma lista de 3 em 3 caracteres
	for i in range (0, len(orf), 3):
		bases.append(orf[i:i+3])
# traduzindo a ORF mais longa em peptídeo, usando a tabela de tradução do início do script
	for base in bases:
		if len(base)==3:
			pep.append(translation_table[base])

# imprimindo a ORF mais comprida em um novo arquivo
# antes, podemos conferir se houve alguma ORF registrada para cada sequência do arquivo
# caso não haja nenhuma, a lista longest_orf continuará com seu valor inicial '0'

	if orf != '0':
		with open ('ORF.fna', 'a') as f:
			print('>'+longest_orf[1]+'\n'+orf, file = f)
	
# imprimindo o peptídeo mais comprido em um novo arquivo
	if orf != '0':
		with open ('ORF.faa', 'a') as g:
			print('>'+longest_orf[1]+'\n'+''.join(pep), file=g)

