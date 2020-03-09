from collections import deque

def test(file1, file2):
    chunkSize = 1000000 # number of characters in each chunk (you will need to adjust this)
    chunk1 = deque([""]) #buffered lines from 1st file
    chunk2 = deque([""]) #buffered lines from 2nd file
    with open(file1, "r") as f1, open(file2, "r") as f2:
        while chunk1 and chunk2:
            line_f1 = chunk1.popleft()
            if not chunk1:
                line_f1,*more = (line_f1+f1.read(chunkSize)).split("\n")
                chunk1.extend(more)
						line_f2 = chunk2.popleft()
            if not chunk2:
                line_f2,*more = (line_f2+f2.read(chunkSize)).split("\n")
                chunk2.extend(more)
            print(line_f1[-1])
            print(line_f2[-1])

if __name__ == "__main__":

    file1 = "/home/rothlab1/rli/dev/tilseq_mutcount/output/test_CBS_170818_2020-02-10-22-35-28/ds_sam_files/128_R1_001_ds.sam"
    file2 = "/home/rothlab1/rli/dev/tilseq_mutcount/output/test_CBS_170818_2020-02-10-22-35-28/ds_sam_files/128_R2_001_ds.sam"
    test(file1, file2)
