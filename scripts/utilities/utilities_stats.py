import math,random,time
from cmath import *
import decimal as dc


from itertools import groupby

g_gamma = 7
p_gamma = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7]
eulers_e = 2.71828182845904523

def gini_coefficient(column):
	diffsum = 0
	try:
		for i, xi in enumerate(column[:-1], 1):
			diffsum += sum([abs(xi -y) for y in column[i:]])

		return  diffsum / (len(column)**2 * sum(column)/len(column))
	except:
		return 0

def mann_whitney_U(sample,control):
	st = time.time()
	#print "#"
	ranked = (control + sample)
	ranked.sort()

	n1= len(control)
	n2 = len(sample)

	r1 = 0
	r2 = 0

	m_u = float(n1*n2)/2

	#print n1, n2
	std_u = n1*n2
	std_u = std_u*(n1+ n2 + 1)
	std_u = float(std_u)/12
	std_u = math.sqrt(std_u)

	rankIndex = {}
	rankCount = dict([(x,len(list(y))) for x,y in groupby(sorted(ranked))])
	ranked.reverse()
	rankIndex = dict([(ranked[x],len(ranked) - x) for x in  range(0,len(ranked))])

	lastTerm = ""

	#print std_u

	control.sort()
	for term in control:
		rankTerm = rankIndex[term]
		countTerm = rankCount[term]


		if countTerm > 1:
			corrected_rank = rankTerm + float(countTerm - 1)/2
		else:
			corrected_rank = rankTerm

		r1 += corrected_rank


	lastTerm = ""
	sample.sort()

	for term in sample:
		rankTerm = rankIndex[term]
		countTerm = rankCount[term]

		if countTerm > 1:
			corrected_rank = rankTerm + float(countTerm - 1)/2
		else:
			corrected_rank = rankTerm

		r2 += corrected_rank

	#print "3",time.time() - st
	T =  tiecorrection(ranked)
	#print T
	#print "4",time.time() - st

	if T == 0:
		print("Error")
		return 1
	else:
		u1 = ((T*n1*n2) + float(n1*(n1+1))/2) - r1
		u2 = ((T*n1*n2) + float(n2*(n2+1))/2) - r2

		z1 = float(u1 - m_u)/std_u
		z2 = float(u2 - m_u)/std_u

		#print u1,u2,z1,z2

		return normal(z1,0,1)

def tiecorrection(sorted):

	n = len(sorted)
	T = 0.0
	i = 0

	while (i<n-1):
		if sorted[i] == sorted[i+1]:
			nties = 1
			while (i<n-1) and (sorted[i] == sorted[i+1]):
				nties = nties +1
				i = i +1
			T = T + nties**3 - nties
		i = i+1
	T = T / float(n**3-n)
	return 1.0 - T

def std_dev(list_temp):
	try:
		mean = sum(list_temp)/len(list_temp)

		sd = 0

		for v in list_temp:
			sd += (v**2 - mean**2)

		sd = sd/len(list_temp)
		return math.sqrt(sd)
	except:
		return -10

def gamma(z):
	z = complex(z)

	if z.real < 0.5:
		return pi / (sin(pi*z)*gamma(1-z))
	else:
		z -= 1
		x = p_gamma[0]
		for i in range(1, g_gamma+2):
			x += p_gamma[i]/(z+i)
		t = z + g_gamma + 0.5

		return sqrt(2*pi) * t**(z.real+0.5) * exp(-t.real) * x

def incomplete_gamma(s,x,alternate=False):
	#Described in computing the incomplete Gamma function to arbitrary precision  - Serge Winitzki
	summer = 0


	if alternate:
		fast = False
		if s > 0:
			for n in range(0,150):
				try:
					if fast:
						#using log - not as accurate
						num = log((x**(n+s))).real
						den = factorialRamanujan(n,logged=True)
						bit = (1.0/(s+n))*((-1)**n)*(e**(num - den))
					else:
						num = ((-1)**n)*(x**(n+s)).real
						den = ((s + n)*factorial(n))
						bit = num/den




					try:
						summer += bit
					except:
						print(("Error", s,x,num,den))
				except:
					pass

		return (gamma(s) - summer).real
	else:
		i = factorial(s-1)
		j = e**(-x)
		l = 0
		for k in range(0,s):
			l +=  float((x**k).real)/factorial(k)

		return  (i*j*l).real

def factorialGamma(m,s=0):
	if m == 0:
		return 1
	else:
		return gamma(m+1).real

def factorialRamanujan(m,s=0,logged=True):
	if m == 0:
		if logged:
			return log(1)
		else:
			return 1


	else:
		if logged:
			return (m*log(m)) - m + (log(m*(1+4*m*(1+2*m)))/6) + (log(pi)/2).real
		else:
			return int(math.ceil(e**((m*log(m)) - m + (log(m*(1+4*m*(1+2*m)))/6) + (log(pi)/2)).real))

def factorial(m,s=0):
	value = int(1)
	if m != 0:
		while m != s:
			value = value*m
			m = m - 1
	return value

def fishers(a,b,c,d):
	x = (factorialRamanujan(a+b) + factorialRamanujan(c+d) + factorialRamanujan(a+c) + factorialRamanujan(b+d))
	y = (factorialRamanujan(a) + factorialRamanujan(b) + factorialRamanujan(c) + factorialRamanujan(d) + factorialRamanujan(a+b+c+d))

	return e**(x - y).real

def poisson(expected,observed):
	expected = int(expected)

	try:
		x = math.exp(-expected)
		y =  pow(expected,observed)

		try:
			z =  math.log(factorial(observed))
		except:
			z =  factorialRamanujan(observed,logged=True).real

		if observed == 0:
			lnExpected = 1
		else:
			lnExpected = math.log(expected)

		return math.exp(observed*lnExpected - expected - z)

	except Exception as e:
		return 0.0


def cum_poisson(expected,observed):
	prob = 1.0
	for x in range(0,observed  +1):
		prob -= poisson(expected,x)


	if prob >= 0:
		return prob
	else:
		return 0.0

def binomial(k,n,p):
	if n > 140:
		return cum_poisson(int(n*p),k)
	else:
		c_bit = factorialRamanujan(n,logged=False)/(factorialRamanujan(k,logged=False)*factorialRamanujan(n-k,logged=False))
		p_bit = (p**k)*((1-p)**(n-k))
		return c_bit * p_bit

def cum_binomial(k,n,p):
	summer = 0
	for i in range(k,n+1):
		summer += binomial(i,n,p)
	return summer

def uniform_product_density(p,n):
	a = (-1)**(n-1)
	b = factorialRamanujan(n-1)
	c = log(p).real**(n-1)

	return (a/b)*c

def cum_uniform_product(p,n,alternate = False):
	try:
		if p == 1:
			return 1
		else:
			if alternate:

				a = (-1)**(n)
				b = incomplete_gamma(n,-log(p))
				c = (log(p).real)**n
				d = factorialRamanujan(n-1,logged=False)
				e = (-log(p).real)**n
				return (a*b*c)/(d*e)

			else:
				a = incomplete_gamma(n,-log(p))
				b = factorial(n-1)
				return a/b
	except Exception as e:
		print(("Error in Sig correction:",p,n))
		print(e)
		return -1

def erf(z):
	summer = 0
	for i in range(0,100):

		try:
			num = ((-1.0)**i) * z**(2*i+1)
		except OverflowError:
			print("Overflow")

		denum = factorial(i)*(2*i+1)
		summer += num/denum


	erf = summer*(2 /math.sqrt(math.pi))
	if erf > 1 or erf < -1:
		erf = 1

	return erf

def choose(k,n):
	return factorial(n)/(factorial(k)*factorial(n-k))

def product(values):
	#prod = 1
	prod = dc.Decimal(str(1))
	for v in values:
		prod *= dc.Decimal(str(v))

	return prod


def normal(z,u,sd):
	if z > 6.2:
		return 0.00000000001
	if z < -6.2:
		return 1
	else:
		return 1 - (1 +erf((z - u)/(abs(sd)*math.sqrt(2))))/2

def rlcProb(rlc,u,sd):
	if rlc > 6.2:
		return 0.00000000001
	if rlc < -6.2:
		return 1
	else:
		return 1 - (1 +erf((rlc - u)/(abs(sd)*math.sqrt(2))))/2

if __name__ == "__main__":
	st = time.time()

	sample = [0.80, 0.83, 1.89, 1.04, 1.45,1.38, 1.91, 1.64, 0.73, 1.46]
	control = [1.15, 0.88, 0.90, 0.74,1.21]

	print((mann_whitney_U(sample,control)))
	from scipy.stats import mannwhitneyu
	U1, p = mannwhitneyu(sample, control,method="exact",alternative='greater')
	print(p )

	sample = [85,89,86,91,77,93,100,82,92,86,86]
	control = [83,73,65,65,90,77,78,97,85,75]

	print((mann_whitney_U(sample,control)))
	from scipy.stats import mannwhitneyu
	U1, p = mannwhitneyu(sample, control,method="exact",alternative='greater')


	print((time.time() - st))
	"""
	# test cum_uniform_product
	summer = 1


	binDict = {}
	samples = 100000
	for i in range(0,samples):
		randInt = int(log(random.random()*random.random()*random.random()*random.random()*random.random(),10).real)
		if randInt in binDict:
			binDict[randInt] += 1
		else:
			binDict[randInt] = 1

	sorter = binDict.keys()
	sorter.sort()
	#sorter.reverse()
	summer = 0

	for bin in sorter:
		summer += binDict[bin]
		print bin,binDict[bin],10**float(bin),float(samples - summer)/samples,1 - cum_uniform_product(10**float(bin),5), cum_uniform_product(10**float(bin),5)
	#"""

	#for i in range(0,10000):
	#	print i,factorial(i)

	#print poisson(2,1.43)
	#print cum_poisson(2,4)
	#print "+++++"
	#print cum_poisson(231,229)
	"""
	for n in range(1,10):
		for  i in range(1,15):
			p = 0.1**i
			print n,"\t",
			print i,"\t",

			print poisson(n,i)
			#print "%1.5g"%(cum_uniform_product(p,n)),"\t",
			#print cum_uniform_product(p,n,alternate=True)

			#print p,-log(p),"%1.3g"%incomplete_gamma(n,-log(p)),"%1.3g"%incomplete_gamma(n,-log(p),True)

		print time.time() - st

		"""
