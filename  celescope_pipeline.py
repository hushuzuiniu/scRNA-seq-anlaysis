# This is a sample Python script.
import os
import argparse
from subprocess import Popen, PIPE, call
import pandas as pd
import logging


# 1:fastq_prefix 2:fastq_dir 3:sample
# fq_r1, fq_r2, sample, fastq_path
def get_fastq(path, all_files):
    file_list = os.listdir(path)
    for file in file_list:
        # 利用os.path.join()
        cur_path = os.path.join(path, file)
        if os.path.isdir(cur_path):
            get_fastq(cur_path, all_files)
        else:
            all_files.append([file, cur_path])
    return all_files


def list_groups(init_list, children_list_len):
    list_of_groups = zip(*(iter(init_list),) * children_list_len)
    end_list = [list(i) for i in list_of_groups]
    count = len(init_list) % children_list_len
    end_list.append(init_list[-count:]) if count != 0 else end_list
    return end_list


def sc_ran_seq(input_dir, output_dir, grch_gt, grh_fa, fastqc_thread, celescope_shell_thread):
    # test_log_level()
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    fastq_list_info = []
    all_fastq_list = get_fastq(input_dir, [])
    for info in all_fastq_list:
        print(info)
        if info[0][-9:] == '.fastq.gz':
            if "_R1" in info[0]:
                fq_r1 = info[0]
                fq_r2 = info[0].replace('_R1', '_R2')
                fastq_r2_path = str(info[1]).replace(info[0], fq_r2)
                if os.path.exists(fastq_r2_path):
                    sample = info[0].split('_R1')[0]
                    fastq_path = str(info[1]).replace(info[0], "")
                    fastq_list_info.append([fq_r1, fq_r2, sample, fastq_path])


    gtf = output_dir + '/gtf'
    # call(, shell=True)
    gtf_command = Popen(["/bin/bash", "-c", 'mkdir -p ' + gtf], stderr=PIPE, stdout=PIPE)
    std_out, std_err = gtf_command.communicate()
    print(std_out)
    print(std_err)
    print("创建gtf完成")
    filtered_gtf = gtf + "/filtered.gtf"

    out_make_ref = output_dir + '/' + 'make_ref'
    # call('mkdir -p ' + out_make_ref, shell=True)
    out_make_ref_command = Popen(["/bin/bash", "-c", 'mkdir -p ' + out_make_ref], stderr=PIPE, stdout=PIPE)
    std_out, std_err = out_make_ref_command.communicate()
    print(std_out)
    print(std_err)
    print("创建make_ref完成")

    map_file = output_dir + "/rna.mapfile"

    multi_rna = output_dir + "/multi_rna"
    # call('mkdir -p ' + multi_rna, shell=True)
    out_make_ref_command = Popen(["/bin/bash", "-c", 'mkdir -p ' + multi_rna], stderr=PIPE, stdout=PIPE)
    std_out, std_err = out_make_ref_command.communicate()
    print(std_out)
    print(std_err)
    print("创建multi_rna完成")

    print("####################celescope utils mkgtf#############################")
    celescope_process_command = "celescope utils mkgtf " + grch_gt + " " + filtered_gtf
    print("celescope_process:" + celescope_process_command)
    celescope_process = Popen(
        ["/bin/bash", "-c", "source activate && conda activate celescope && " + celescope_process_command], stderr=PIPE,
        stdout=PIPE)
    std_out, std_err = celescope_process.communicate()
    print(std_out)
    print(std_err)
    print("##########################celescope utils mkgtf 结束######################")

    print("############################celescope rna mkref 开始#########################")
    # Create a reference genome file with the mkref comman
    mkref_process_command = "celescope rna mkref --genome_name Homo_filtered --fasta " + grh_fa + " --gtf " + filtered_gtf
    print("mkref_process_command:" + mkref_process_command)
    mkref_process = Popen(["/bin/bash", "-c",
                           "cd " + out_make_ref + " && source activate && conda activate celescope && " + mkref_process_command],
                          stderr=PIPE, stdout=PIPE)
    std_out, std_err = mkref_process.communicate()
    print(std_out)
    print(std_err)
    print("##############################celescope rna mkref 结束#########################")

    fastq_list_groups = list_groups(fastq_list_info, fastqc_thread)

    for group in fastq_list_groups:
        print(group)
        print("do fastqc")

        fastqc_command = []
        for info in group:
            fastqc_out_dir = output_dir + '/' + info[2]
            call('mkdir -p ' + fastqc_out_dir, shell=True)
            fastqc_command.append(
                'fastqc ' + info[3] + "/" + info[0] + ' ' + info[3] + "/" + info[1] + ' -o ' + fastqc_out_dir)

        print("################################fastqc 开始##################################")
        print(" & ".join(fastqc_command))
        fastqc_process = Popen(["/bin/bash", "-c", " & ".join(fastqc_command)],
                               stderr=PIPE, stdout=PIPE)
        std_out, std_err = fastqc_process.communicate()
        print(std_out)
        print(std_err)
        print("############################fastqc 结束######################################")

    print(fastq_list_info)
    map_file_list = []
    for fastq in fastq_list_info:
        # 1:fastq_prefix 2:fastq_dir 3:sample
        print(fastq)
        map_file_list.append([fastq[2], fastq[3], fastq[2]])

    df_map_file = pd.DataFrame(map_file_list)
    df_map_file.to_csv(map_file, sep='\t', header=None, index=None)
    print("############################multi_rna 开始######################################")
    multi_rna_command = "source activate && conda activate celescope && cd " + out_make_ref + " && multi_rna --mapfile " + map_file + " --genomeDir " + out_make_ref + " --thread 12 --mod shell --outdir " + multi_rna
    print(multi_rna_command)
    mkref_process = Popen(["/bin/bash", "-c", multi_rna_command],
                          stderr=PIPE, stdout=PIPE)
    std_out, std_err = mkref_process.communicate()
    print(std_out)
    print(std_err)
    print("############################multi_rna 结束######################################")
    for index, row in df_map_file.iterrows():
        print(row) 
        mk_ref_process_chmod_command = out_make_ref + "/shell/" + row[2] + ".sh"

        print("修改执行权限")
        chmod = Popen(["/bin/bash", "-c", "chmod 777 " + mk_ref_process_chmod_command], stderr=PIPE, stdout=PIPE)
        std_out, std_err = chmod.communicate()
        print(std_out)
        print(std_err)
        print("修改权限结束")

    df_map_file.columns = ["A", "B", "C"]
    list_of_map_file = list(df_map_file['C'])
    list_of_map_file_groups = list_groups(list_of_map_file, celescope_shell_thread)

    for rows in list_of_map_file_groups:
        print(rows)
        sample_command = []
        for sample in rows:
            sample_command.append(out_make_ref + "/shell/" + sample + ".sh")
        print("mk_ref_process_command:" + " & ".join(sample_command))
        mk_ref_process = Popen(
            ["/bin/bash", "-c", "source activate && conda activate celescope && " + " & ".join(sample_command)],
            stderr=PIPE, stdout=PIPE)
        std_out, std_err = mk_ref_process.communicate()
        print(std_out)
        print(std_err)


# Define input folder and output folder
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fastq_dir', type=str, default=None)
    parser.add_argument('--output_dir', type=str, default=None)
    parser.add_argument('--grch_gt', type=str, default=None)
    parser.add_argument('--grh_fa', type=str, default=None)
    parser.add_argument('--fastqc_thread', type=int, default=1)  
    parser.add_argument('--celescope_shell_thread', type=int, default=1)  
    args = parser.parse_args()
    sc_ran_seq(args.fastq_dir, args.output_dir, args.grch_gt, args.grh_fa, args.fastqc_thread,
               args.celescope_shell_thread)
