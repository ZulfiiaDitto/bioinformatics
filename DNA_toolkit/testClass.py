from bioseq import bio_seq
from utilities import writeTextFile, readTextFile

# testing for the  class structure

testDna = bio_seq()
#print(testDna.show_seq_info())
testDna.generate_rand_seq(length=40,seq_type= "RNA")
print(testDna.countNucFrequency())
print(testDna.transcription())
print(testDna.gc_content())
print(testDna.gc_content_subset())
print(testDna.translate_seq())
print(testDna.codon_usage('T'))
for i in testDna.gen_reading_frames():
    print(i)

print( testDna.proteins_from_rf(['L', 'R', 'V', 'M', 'A', 'Y', 'C', 'M', 'L', 'E', 'Y', '_', 'G']))
print( testDna.all_proteins_from_all_pr_frames(end = 30))

writeTextFile('test.txt', testDna.seq)