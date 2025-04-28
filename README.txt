For general documentation see README_v1.

See WHATSNEW for latest information

Version 2.0

* has Ilink capability
* handles genotype specific penetrances.
  Specify 4 penetrances instead of 3. The order of the trait genotypes
  is  1|1   1|2   2|1   2|2  with paternal allele on the left. Note
  that the order of the heterozygotes is opposite that of ghi (genehunter
  imprinting version). To run imprinting add the command line option
  -L I.
* uses the child vector algorithm as described in the reference below.
  VITESSE will automatically switch between the parental pair and child
  vector algorithm. To run only the parental pair algorithm (original
  vitesse) use the command line option -1 f. To run only the child vector
  algorithm use the command line option -1 c.
* does NOT do the byfamily likelihoods (lods) by default as version 1 
  does. To obtain the byfamily results use the command line option -F.
* can print the lod scores for a LINKMAP run by using the command
  line option -M
* can use non-default names for pedigree and data files by using the
  command line options -p <pedfile name> and -d <datafile name>, respectively.
* accepts phased data in the pedigree file. Genotypes in LINKAGE format
  are written as A B, which a priori represents 2 phases, A | B and B | A
  where A | B represents the genotype with A being the paternal and B the 
  maternal allele. VITESSE accepts genotypes of the form A | B.
* accepts groups of genotypes at a locus. Enclose the group with { }. Example
  { 1 2  3 4 } represents two genotypes 1 2 and 3 4. Thus { A | B  B | A } 
  is the same as A B.
* accepts half-typed genotypes. If one allele is known with certainty, say 
  the A allele, then A 0 can be entered in the pedigree file. Normally VITESSE
  treats this as an error. To override the error message use the command line 
  option -Z.


Complexity:
Certain data structures in VITESSE are optimized for up to the 10 marker range,
since that was the range that was of interest for large pedigrees with missing
data. Beyond 10 markers these data structures become very large.
Thus the runtime complexity above 10 markers does not always reflect the
intrinsic complexity of a problem. For example, a fully-typed pedigree with
20 markers may be less complex than a pedigree with 3 markers and mostly missing
data. VITESSE is optimized for the latter, but not the former.

If you use this version of VITESSE please also cite:

O'Connell JR (2001) Rapid Multipoint Linkage Analysis via Inheritance Vectors 
in the Elston-Stewart Algorithm.  Hum Hered. 51(4):226-40.

