'''
This file parses the log files such as build__2_1.log
'''

from collections import defaultdict
import re

def create_summary(file):
    '''
        This function parses the test results from the build__2_1.log file
        into a dictionary
    '''
    sum_data={}
    with open(file,"r") as fp:
        lines=fp.read().splitlines()
        # Find the line where the summary starts
        i=0
        while not re.match("^\s*Summary for configuration sim",lines[i]):
            i+=1
        i+=1
        # Loop until the end of the summary
        while not re.match("\s*Tests passed:",lines[i]) and not re.match("\s*Tests failed:",lines[i]):
            # The spacing of this line is unique and as such requires a special if statement
            if re.match("\s*Number passed only to\s*",lines[i]) and lines[i]!="":
                split_l=lines[i+1].split("->")
                split_l[0]=lines[i]+" "+split_l[0]

                # Convert numerical data to integer
                try:
                    sum_data[" ".join(split_l[0].split())]=int(split_l[1].strip())
                except:
                    sum_data[" ".join(split_l[0].split())]=split_l[1].strip()
                # This data field has to lines and as such increment the line counter twice
                i+=1

            # The regular line parsing code
            elif lines[i]!="":
                split_l=lines[i].split("->")
                try:
                    val = int(split_l[1])
                except ValueError:
                    val = split_l[1].strip()
                sum_data[" ".join(split_l[0].split())]=val
            i+=1
            
    return sum_data

def get_tests(readfile):
    '''
        This functions reads the names of the tests that failed and passed
        from the log files and returns the set of them in tuples.
    '''
    passed=set()
    failed=set()
    with open(readfile,"r") as fp:
        line = fp.readline()
        while line:
            if re.match("\s*Tests passed:",line):
                line = fp.readline()
                while line:
                    line = fp.readline()
                    if line.strip() == "":
                        break
                    passed.add(" ".join(line.split()))
            elif re.match("\s*Tests failed:",line): 
                line = fp.readline()
                while line:
                    line = fp.readline()
                    if line.strip() == "":
                        break
                    failed.add(" ".join(line.split()))
            line = fp.readline()
    return passed,failed

def test_comp(readfile_new,readfile_old):
    '''
        This function compares the tests from the previous run 
        to see which failed, which failed this time but did
        not in the last, which passed this time but not in the last,
        and the new and removed tests.
    '''
    passed_n,failed_n=get_tests(readfile_new)
    try:
        passed_o,failed_o=get_tests(readfile_old)
    except FileNotFoundError:
        passed_o,failed_o=(set(), set())
    newly_p=passed_n-passed_o
    newly_f=failed_n-failed_o
    new_tests=(passed_n.union(failed_n))-(passed_o.union(failed_o))
    missing_tests=(passed_o.union(failed_o))-(passed_n.union(failed_n))
    def test_key(label):
        # Cactus test names allow all characters acceptable in file names,
        # though ones with spaces would be rare
        # Cactus thorn names must be C identifiers
        m = re.match("([^ ]*) [(]from (\w*)[)]", label)
        if m:
            return (m.group(2), m.group(1))
        else:
            return ('', label)
    types=["Failed Tests","Newly Passing Tests","Newly Failing Tests","Newly Added Tests", "Removed Tests"]
    tests=[failed_n,newly_p,newly_f,new_tests,missing_tests]
    sorted_tests=[sorted(test, key=test_key) for test in tests]
    test_dict={types[i]:sorted_tests[i] for i in range(len(types))}
    return test_dict

def get_times(readfile):
    '''
        This function finds the times taken for each test in the log
        file and then stores that in a dictionary and then sorts
        those tests in descending order by time
    '''
    times={}
    with open(readfile,"r") as fp:
        lines=fp.read().splitlines()
        ind=0
        for line in lines:
            if re.match("\s*Details:",line):
                break
            ind+=1
        while ind<len(lines):
            try:
                time_i=lines[ind].index('(')
                tim=float(lines[ind][time_i+1:].split()[0])
                test_name=lines[ind][:time_i-1].split()[0]
                times[test_name]=tim
            except:
                pass
            ind+=1

    return {test:ti for test,ti in sorted(times.items(),key= lambda x : x[1],reverse=True)} # This is a dictionary comprehension that uses sorted to order the items in times.items() into a dictionary
def exceed_thresh(time_dict,thresh):
    '''
        This function finds tests that exceed a certain time threshhold
    '''
    return {test:ti for test,ti in time_dict.items() if ti>thresh}

def longest_tests(time_dict,num_tests):
    '''
        This function finds the tests that took the longest time.
        It returns a dictionary of the num_test number of longest
        tests.
    '''
    longest={}
    i=0
    for item in time_dict.items():
        if i>num_tests-1:
            break
        longest[item[0]]=item[1]
        i+=1
    return longest

def get_unrunnable(readfile):
    '''
        This test reads the logfile looking for tests that could not be run
        and the corresponding reason.
    '''
    miss_th={}
    miss_proc={}
    with open(readfile,"r") as fp:
        line = fp.readline()
        while line:
            if re.match("\s*Tests missed for lack of thorns:",line):
                line = fp.readline() # read empty line
                while line:
                    line = fp.readline()
                    if not re.match("\s*\w+ in \w+", line):
                        break
                    thorn = line.split()[0]
                    fp.readline()
                    line = fp.readline()
                    missing = line.split(":")[1].split()
                    miss_th[thorn] = missing
                line = fp.readline()
                continue # re-parse line
            elif re.match("\s*Tests missed for different number of processors required:",line): 
                line = fp.readline() # read empty line
                while line:
                    line = fp.readline()
                    if not re.match("\s*\w+ in \w+", line):
                        break
                    thorn = line.split()[0]
                    fp.readline()
                    line = fp.readline()
                    missing = line.split(":")[1].split()
                    miss_proc[thorn] = missing
                    line = fp.readline()
                continue # re-parse line
            line = fp.readline()
    return miss_th,miss_proc

def get_data(name):
    '''
        Retrieves singular field of data from the data csv as a list
    '''
    data={}
    with open('test_nums.csv','r') as csvfile:
        fields=csvfile.readline().strip().split(",")
        name_i=fields.index(name)
        line=csvfile.readline()
        while line !="":
            entry=line.strip().split(",")
            data[entry[0]]=float(entry[name_i])
            line=csvfile.readline()
    return data

def get_compile(name):
    '''
        Gets the total number of compile time warnings by using
        the get_warning_thorns function
    '''
    return sum(get_warning_thorns(name).values())

def get_warning_type(name):
    '''
        Compiles of counts of what types of warnings are produced the most
        during compilation
    '''
    warning_types=defaultdict(int)
    
    with open(name) as build:
        for line in build.readlines():
            m = re.search(".*/sim/build/([^/]*).* [wW]arning:", line)
            if(m):
                warning_types[line[line.find("[-W"):-1]]+=1
    return warning_types

def get_warning_thorns(name):
    '''
        This code finds how many compile time warnings are related each thorn
    '''
    warning_types=defaultdict(int)
    i=0
    count=0
    with open(name) as build:
        lines=build.readlines()
        for line in lines:
            i+=1
            # This regex search finds inline warnings based on the pattern given
            inline = re.search(".*/sim/build/([^/]*).* [wW]arning:", line)

            # This regex search finds the pattern shown below as twoline warnings are structure in this way
            twoline= re.search("[wW]arning:.*at",line)

            if(inline):

                trunc=line[line.find("build/")+6:-1]
                trunc=trunc[:trunc.find("/")]
                warning_types[trunc]+=1
            if(twoline):
                count+=1
                nextl=lines[i+1]
                nextnextl=lines[i+2]
                warning=re.search(".*/sim/build/([^/]*).*", nextl)
                warning2=re.search(".*/sim/build/([^/]*).*", nextnextl)
                if(warning):
                    trunc=nextl[nextl.find("build/")+6:-1]
                    trunc=trunc[:trunc.find("/")]
                    warning_types[trunc]+=1
                if(warning2):
                    trunc=nextnextl[nextnextl.find("build/")+6:-1]
                    trunc=trunc[:trunc.find("/")]
                    warning_types[trunc]+=1
    return warning_types
