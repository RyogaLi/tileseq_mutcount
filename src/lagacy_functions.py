## functions that are no longer in use
import os
import random

def downsample(n, r1, r2, output_path):
    """
    Randomly downsample r1 and r2 to n reads
    """

    record_number = 0

    base_r1 = os.path.basename(r1)
    base_r2 = os.path.basename(r2)

    output_r1 = output_path+"/"+base_r1.replace(".fastq", "_ds.fastq")
    output_r2 = output_path+"/"+base_r2.replace(".fastq", "_ds.fastq")

    with open(r1, "r") as i:
        r1_lines = sum([1 for line in i])
    with open(r2, "r") as j:
        r2_lines = sum([1 for line in j])

    total_records_r1 = int(r1_lines/4)
    total_records_r2 = int(r2_lines/4)

    if n > total_records_r1 or n > total_records_r2:
        os.system("cp "+r1+ " "+output_r1)
        os.system("cp "+r2+ " "+output_r2)
        return output_r1, output_r2
    else:
        if total_records_r1 <= total_records_r2:
            records_to_keep = set(random.sample(range(total_records_r1 + 1), n))
        else:
            records_to_keep = set(random.sample(range(total_records_r2 + 1), n))

    with open(r1, "r") as i_1:
        with open(output_r1, "w") as o1:

            for line1 in i_1:
                line2 = i_1.readline()
                line3 = i_1.readline()
                line4 = i_1.readline()
                if record_number in records_to_keep:
                    o1.write(line1)
                    o1.write(line2)
                    o1.write(line3)
                    o1.write(line4)
                record_number += 1


    record_number = 0
    with open(r2, "r") as i_2:
        with open(output_r2, "w") as o2:
            for line1 in i_2:
                line2 = i_2.readline()
                line3 = i_2.readline()
                line4 = i_2.readline()
                if record_number in records_to_keep:
                    o2.write(line1)
                    o2.write(line2)
                    o2.write(line3)
                    o2.write(line4)
                record_number += 1

    return output_r1, output_r2
