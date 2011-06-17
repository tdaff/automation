#! /usr/bin/env python

import os,sys,subprocess

def exit():
   print "Exiting."
   sys.exit()

def trim_non_dump_lines(file):
    temp_list = file.readlines()
    for i in range(len(temp_list)):
        if temp_list[i].find("#dump#") != 0:
            temp_list[i] = ""
    while True:
        try:
            temp_list.remove("")
        except:
            break
    return temp_list

def write_list_to_file(list,filename):
    f1 = open(filename,"w")
    for i in range(len(list)):
        f1.write(list[i])
    f1.close()   

def compare_dumps(dumpfile):
    template_file = open("./Test_Pack/TEMPLATES/"+dumpfile,'r')
    dump_file = open("./dump",'r')
    template = trim_non_dump_lines(template_file)
    dump = trim_non_dump_lines(dump_file)
    dump_file.close()
    template_file.close()
    write_list_to_file(template,"template_tmp")
    write_list_to_file(dump,"dump_tmp")
    lines_different = subprocess.popen("diff dump_tmp template_tmp| wc -l",stdout=subprocss.PIPE,Shell=True).stdout()
    if lines_different == 0:
       print "Everything is being generated perfectly.\n Test passed."

def simple_test():
    os.system("cp ./Test_Pack/INPUTS/JOB_INPUT.1 ./JOB_INPUT")
    os.system("00MN_AUTO_OCC INPUT_MOF1 mof1systest > dump")
    compare_dumps("dump.1")

def regular_test():
    os.system("cp ./Test_Pack/INPUTS/JOB_INPUT.2 ./JOB_INPUT")
    os.system("00MN_AUTO_OCC INPUT_MOF1 mof1systest > dump")
    compare_dumps("dump.2")

def heavy_test():
    os.system("cp ./Test_Pack/INPUTS/JOB_INPUT.3 ./JOB_INPUT")
    os.system("00MN_AUTO_OCC INPUT_MOF1 mof1systest > dump")
    compare_dumps("dump.3")

def extreme_test():
    os.system("cp ./Test_Pack/INPUTS/JOB_INPUT.4 ./JOB_INPUT")
    os.system("00MN_AUTO_OCC INPUT_MOF1 mof1systest > dump")
    compare_dumps("dump.4")


def main(argv):
    OPTIONS = {1:simple_test,2:regular_test,3:heavy_test,4:extreme_test,5:exit}

    menu = "-----------------\n(1) Perform a Simple System Test\n(2) Perform a Regular System Test\n(3) Perform a Heavy System Test\n(4) Perform an EXTREME System Test\n(5) Exit\n-----------------"

    print menu
    choice = raw_input("Enter an Option: ")

    valid_list = ["1","2","3","4","5"]

    while True:
        choice_valid = False
        while not(choice_valid):
            if valid_list.__contains__(choice):
                choice_valid = True
            else:
                print "Invalid option: please choose again\n"
                print menu
                choice = raw_input("Enter an Option:")
        OPTIONS[int(choice)]()
        print menu
        choice = raw_input("Enter an Option:")
    




if __name__ == "__main__":
    sys.exit(main(sys.argv))

